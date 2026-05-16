use std::collections::BTreeSet;

include!(concat!(env!("OUT_DIR"), "/hll_bias.rs"));

/// MurmurHash3 finalizer variant used by HyperLogLog.
/// Note: adds 1 to key before hashing (to avoid hash(0) == 0).
#[inline]
fn murmurhash3_finalizer(mut key: u64) -> u64 {
    key = key.wrapping_add(1);
    key ^= key >> 33;
    key = key.wrapping_mul(0xff51afd7ed558ccd);
    key ^= key >> 33;
    key = key.wrapping_mul(0xc4ceb9fe1a85ec53);
    key ^= key >> 33;
    key
}

/// Manual count-leading-zeros for u64 using the Hacker's Delight binary
/// search algorithm. Kept as a portable fallback; production code should
/// prefer `u64::leading_zeros` (compiled to `lzcnt`/`clz`).
#[allow(dead_code)]
fn clz_manual(mut x: u64) -> i32 {
    let mut n = 64;
    let mut y = x >> 32;
    if y != 0 {
        n -= 32;
        x = y;
    }
    y = x >> 16;
    if y != 0 {
        n -= 16;
        x = y;
    }
    y = x >> 8;
    if y != 0 {
        n -= 8;
        x = y;
    }
    y = x >> 4;
    if y != 0 {
        n -= 4;
        x = y;
    }
    y = x >> 2;
    if y != 0 {
        n -= 2;
        x = y;
    }
    y = x >> 1;
    if y != 0 {
        return n - 2;
    }
    n - x as i32
}

/// Count leading zeros of a `u32`, returning `max` when `x == 0`
/// (matching the C++ `clz(x, max)` helper that wraps `__builtin_clz`).
#[inline]
fn clz_u32(x: u32, max: u8) -> u8 {
    if x == 0 {
        max
    } else {
        x.leading_zeros() as u8
    }
}

/// Count leading zeros of a `u64`, returning `max` when `x == 0`
/// (matching the C++ `clz(x, max)` helper that wraps `__builtin_clzl`).
#[inline]
fn clz_u64(x: u64, max: u8) -> u8 {
    if x == 0 {
        max
    } else {
        x.leading_zeros() as u8
    }
}

/// Branchless variant of `clz` that injects a sentinel `1` bit before the
/// shift so the result is bounded. Mirrors C++ `clz_p` (unused alternative).
#[allow(dead_code)]
#[inline]
fn clz_p(x: u32, p: u8) -> u8 {
    (((x << 1) | 1) << (p - 1)).leading_zeros() as u8
}

/// Extract a contiguous bitfield from `value` (LSB-0 numbering, bits `lo..hi`
/// inclusive on the low side, exclusive on the high side). If `shift_left`
/// is true the field is left-aligned to the top of the word; otherwise it
/// is right-aligned.
#[allow(dead_code)]
#[inline]
fn extract_bits(value: u64, hi: u8, lo: u8, shift_left: bool) -> u64 {
    let bitmask = ((1u64 << (hi - lo)) - 1) << lo;
    let result = value & bitmask;
    if shift_left {
        result << (64 - hi)
    } else {
        result >> lo
    }
}

/// Return the top `hi` bits of a 64-bit value, right-justified.
#[allow(dead_code)]
#[inline]
fn extract_high_bits_u64(bits: u64, hi: u8) -> u64 {
    bits >> (64 - hi)
}

/// Return the top `hi` bits of a 32-bit value, right-justified.
#[allow(dead_code)]
#[inline]
fn extract_high_bits_u32(bits: u32, hi: u8) -> u32 {
    bits >> (32 - hi)
}

/// Take the first `p` bits of a 64-bit hash and use them as the register
/// index. Matches `getIndex(hash_value, p)` from the C++ implementation.
#[inline]
fn get_index_u64(hash_value: u64, p: u8) -> u32 {
    (hash_value >> (64 - p)) as u32
}

/// Take the first `p` bits of a 32-bit hash (sparse-mode encoded value) and
/// use them as the register index.
#[allow(dead_code)]
#[inline]
fn get_index_u32(hash_value: u32, p: u8) -> u32 {
    hash_value >> (32 - p)
}

/// Build a `u64` whose lowest `p` bits are 1.
#[allow(dead_code)]
#[inline]
fn trailing_ones(p: u8) -> u64 {
    (1u64 << p) - 1
}

/// Compute the HLL rank of a 64-bit hash with precision `p`: shift off the
/// `p` index bits, then count leading zeros + 1 over the remaining bits.
#[inline]
fn get_rank_u64(hash_value: u64, p: u8) -> u8 {
    let rank_bits = hash_value << p;
    clz_u64(rank_bits, 64 - p) + 1
}

/// Decode the rank from a sparse-mode 32-bit encoded hash. When the
/// low-bit flag is set, bits between `p` and `p_prime` were known to be
/// zero and the rank stored in bits 1..7 is the additional rank
/// `p_prime - p` above that baseline; otherwise the rank is computed
/// directly from the encoded value's high bits.
#[inline]
fn get_encoded_rank(encoded_hash_value: u32, p_prime: u8, p: u8) -> u8 {
    if (encoded_hash_value & 1) == 1 {
        let additional_rank = p_prime - p;
        let extracted = ((encoded_hash_value >> 1) & 0x3F) as u8; // 6 bits
        additional_rank + extracted
    } else {
        get_rank_u32(encoded_hash_value, p)
    }
}

/// 32-bit variant of `get_rank_u64`, used when decoding from sparse mode.
#[inline]
fn get_rank_u32(hash_value: u32, p: u8) -> u8 {
    let rank_bits = hash_value << p;
    clz_u32(rank_bits, 32 - p) + 1
}

/// Encode a 64-bit hash as a 32-bit sparse-mode value (see Heule et al.
/// 2013, section 5.3). The top `p_prime` bits are placed at the top of
/// the 32-bit result. If bits `p..p_prime` are all zero we lose information
/// about the rank by truncating, so we set the low flag bit and pack the
/// extra rank (`pPrime - p`) into bits 1..7.
#[inline]
fn encode_hash_in_32bit(hash_value: u64, p_prime: u8, p: u8) -> u32 {
    let idx = (extract_high_bits_u64(hash_value, p_prime) as u32) << (32 - p_prime);
    if (idx << p) == 0 {
        let additional_rank = get_rank_u64(hash_value, p_prime);
        idx | ((additional_rank as u32) << 1) | 1
    } else {
        idx
    }
}

/// Insert an encoded hash into a sorted-vec sparse list, keeping only the
/// most informative encoding per index. When two encodings share an index,
/// prefer the one with the flag bit set (more accurate rank); otherwise
/// keep the smaller or larger encoding depending on the flag. Mirrors the
/// C++ `addHashToSparseList(vector&, ...)` overload (currently unused — we
/// use a `BTreeSet`-backed sparse list instead, but the variant is kept
/// for parity with the original implementation).
#[allow(dead_code)]
fn add_hash_to_sparse_vec(vec: &mut Vec<u32>, val: u32, p_prime: u8) {
    let pos = match vec.binary_search(&val) {
        Ok(_) => return,
        Err(pos) => pos,
    };
    if pos == vec.len() {
        vec.insert(pos, val);
        return;
    }

    let current = vec[pos];
    if extract_high_bits_u32(val, p_prime) == extract_high_bits_u32(current, p_prime) {
        if (current & 1) == (val & 1) {
            if (val & 1) == 1 {
                if val > current {
                    vec[pos] = val;
                }
            } else if val < current {
                vec[pos] = val;
            }
        } else if (val & 1) == 1 {
            vec[pos] = val;
        }
    } else {
        vec.insert(pos, val);
    }
}

/// Insert an encoded hash into the BTreeSet-backed sparse list. The set
/// transparently handles de-duplication; we do not currently merge
/// entries that share an index but differ in precision (matches the
/// C++ `addHashToSparseList(SET&, ...)` template specialisation).
#[inline]
fn add_hash_to_sparse_set(set: &mut BTreeSet<u32>, val: u32, _p_prime: u8) {
    set.insert(val);
}

/// HyperLogLog bias-correction constant `alpha_m` as a function of the
/// number of registers. The special cases for 16/32/64 come from the
/// original Flajolet et al. paper; for `m >= 128` the formula
/// `0.7213 / (1 + 1.079/m)` is used.
#[allow(dead_code)]
fn alpha(m: u32) -> f64 {
    match m {
        16 => 0.673,
        32 => 0.697,
        64 => 0.709,
        _ => 0.7213 / (1.0 + 1.079 / m as f64),
    }
}

/// Linear-counting estimator of Whang et al. (1990):
/// `n_hat = m * ln(m / v)` where `m` is the number of bins and `v` is the
/// number of empty bins. Used at low cardinalities where the raw
/// HyperLogLog estimate has large relative error.
#[allow(dead_code)]
fn linear_counting(m: u32, v: u32) -> f64 {
    assert!(v <= m, "number of v should not be greater than m");
    (m as f64) * ((m as f64) / (v as f64)).ln()
}

/// Raw HyperLogLog estimate: `alpha_m * m^2 / sum(2^-M[i])` — the harmonic
/// mean of the register values, scaled by `alpha_m * m`.
#[allow(dead_code)]
fn calculate_raw_estimate(registers: &[u8]) -> f64 {
    let m = registers.len();
    let mut inverse_sum: f64 = 0.0;
    for &r in registers {
        inverse_sum += 1.0 / (1u64 << r) as f64;
    }
    alpha(m as u32) * (m * m) as f64 / inverse_sum
}

/// Count how many registers are still zero (used to trigger linear-counting
/// fallback at low cardinality).
#[allow(dead_code)]
fn count_zeros(registers: &[u8]) -> u32 {
    registers.iter().filter(|&&v| v == 0).count() as u32
}

/// Heule et al. (HLL++) empirical bias correction. Looks up the raw estimate
/// in the pre-computed per-precision tables (`RAW_ESTIMATE_DATA` and
/// `BIAS_DATA`, generated from the published reference data) and returns a
/// linearly interpolated bias to subtract from the raw estimate. Out-of-range
/// estimates clamp to the table endpoints.
fn get_estimate_bias(estimate: f64, p: u8) -> f64 {
    let raw_estimate_table = RAW_ESTIMATE_DATA[(p - 4) as usize];
    let bias_table = BIAS_DATA[(p - 4) as usize];

    if raw_estimate_table[0] >= estimate {
        return bias_table[0];
    }
    if raw_estimate_table[raw_estimate_table.len() - 1] <= estimate {
        return bias_table[bias_table.len() - 1];
    }

    let pos = raw_estimate_table.partition_point(|&value| value < estimate);
    let e1 = raw_estimate_table[pos - 1];
    let e2 = raw_estimate_table[pos];
    let c = (estimate - e1) / (e2 - e1);
    bias_table[pos - 1] * (1.0 - c) + bias_table[pos] * c
}

/// Build the "register histogram" `C` used by Ertl's improved estimator:
/// `C[i]` is the number of registers equal to value `i`. The histogram has
/// size `q + 2` where `q = 64 - p` is the maximum possible rank value.
fn register_histogram(registers: &[u8], q: usize) -> Vec<i32> {
    let mut c = vec![0i32; q + 2];
    for &r in registers {
        c[r as usize] += 1;
    }
    c
}

/// Register histogram constructed from the sparse-mode list rather than the
/// dense register array. All not-yet-observed registers contribute to
/// `C[0]` (the count of zero-valued registers).
fn sparse_register_histogram(
    sparse_list: &BTreeSet<u32>,
    p_prime: u8,
    p: u8,
    q: usize,
) -> Vec<i32> {
    let mut c = vec![0i32; q + 2];
    let mut m = 1usize << p_prime;
    for &encoded_hash_value in sparse_list {
        let rank_val = get_encoded_rank(encoded_hash_value, p_prime, p) as usize;
        c[rank_val] += 1;
        m -= 1;
    }
    c[0] = m as i32;
    c
}

/// Ertl's sigma correction for 0-registers.
fn sigma(mut x: f64) -> f64 {
    assert!((0.0..=1.0).contains(&x));
    if x == 1.0 {
        return f64::INFINITY;
    }
    let mut sigma_x = x;
    let mut y = 1.0;
    loop {
        let prev = sigma_x;
        x *= x;
        sigma_x += x * y;
        y += y;
        if sigma_x == prev {
            break;
        }
    }
    sigma_x
}

/// Alternative formulation of `sigma` from the C++ reference (unused).
/// Terminates when the squared term drops below `f64::EPSILON` rather than
/// when the running sum stops changing.
#[allow(dead_code)]
fn sigma_mod(x: f64) -> f64 {
    assert!((0.0..=1.0).contains(&x));
    if x == 1.0 {
        return f64::INFINITY;
    }

    let mut sigma_x = x;
    let mut x_sq = x * x;
    let mut two_exp = 1.0;
    while x_sq > f64::EPSILON {
        sigma_x += x_sq * two_exp;
        x_sq *= x_sq;
        two_exp += two_exp;
    }
    sigma_x
}

/// Ertl's tau correction for high-value registers.
fn tau(mut x: f64) -> f64 {
    assert!((0.0..=1.0).contains(&x));
    if x == 0.0 || x == 1.0 {
        return 0.0;
    }
    let mut y = 1.0;
    let mut tau_x = 1.0 - x;
    loop {
        let prev = tau_x;
        x = x.sqrt();
        y /= 2.0;
        tau_x -= (1.0 - x).powi(2) * y;
        if tau_x == prev {
            break;
        }
    }
    tau_x / 3.0
}

/// Numerical Recipes random hash. Provided for parity with the C++ source;
/// unused — the cardinality sketch uses `murmurhash3_finalizer`.
#[allow(dead_code)]
fn ranhash(u: u64) -> u64 {
    let mut v = u
        .wrapping_mul(3935559000370003845)
        .wrapping_add(2691343689449507681);
    v ^= v >> 21;
    v ^= v << 37;
    v ^= v >> 4;
    v = v.wrapping_mul(4768777513237032717);
    v ^= v << 20;
    v ^= v >> 41;
    v ^= v << 5;
    v
}

/// Thomas Wang's 64-bit integer mixer. Provided as an alternative bit
/// mixer for HLL hashing; not currently wired in.
#[allow(dead_code)]
fn wang_mixer(mut key: u64) -> u64 {
    key = (!key).wrapping_add(key << 21);
    key ^= key >> 24;
    key = key.wrapping_add(key << 3).wrapping_add(key << 8);
    key ^= key >> 14;
    key = key.wrapping_add(key << 2).wrapping_add(key << 4);
    key ^= key >> 28;
    key = key.wrapping_add(key << 31);
    key
}

/// Identity hash (no mixing). Provided for parity with the C++ `NoHash`
/// functor used in benchmarks; unused.
#[allow(dead_code)]
fn no_hash(u: u64) -> usize {
    u as usize
}

/// HyperLogLogPlusMinus — cardinality estimation for distinct k-mer counting.
/// Port of the C++ HyperLogLogPlusMinus<uint64_t> template.
#[derive(Clone)]
pub struct HyperLogLogPlusMinus {
    p: u8,
    m: usize,
    registers: Vec<u8>,
    n_observed: u64,
    sparse: bool,
    sparse_list: BTreeSet<u32>,
    pub use_n_observed: bool,
}

const P_PRIME: u8 = 25;
const M_PRIME: u32 = 1 << P_PRIME;

impl HyperLogLogPlusMinus {
    /// Construct a sketch with the given precision (number of register bits;
    /// register count `m = 2^precision`). Sparse mode is enabled by default,
    /// matching the C++ constructor default.
    pub fn new(precision: u8) -> Self {
        Self::with_sparse(precision, true)
    }

    /// Construct a sketch with explicit sparse-mode selection. `precision`
    /// must be in `[4, 18]`. When `sparse` is true, the sketch starts in
    /// sparse mode and switches to dense once the sparse list grows past
    /// `m/4` entries.
    pub fn with_sparse(precision: u8, sparse: bool) -> Self {
        assert!(
            (4..=18).contains(&precision),
            "precision must be between 4 and 18"
        );
        let m = 1usize << precision;
        HyperLogLogPlusMinus {
            p: precision,
            m,
            registers: if sparse { Vec::new() } else { vec![0u8; m] },
            n_observed: 0,
            sparse,
            sparse_list: BTreeSet::new(),
            use_n_observed: true,
        }
    }

    /// Deep-copy constructor — Rust analogue of the C++ copy ctor. Prefer
    /// `Clone` in idiomatic Rust code; this is provided for parity.
    pub fn copy_construct(other: &HyperLogLogPlusMinus) -> Self {
        HyperLogLogPlusMinus {
            p: other.p,
            m: other.m,
            registers: other.registers.clone(),
            n_observed: other.n_observed,
            sparse: other.sparse,
            sparse_list: other.sparse_list.clone(),
            use_n_observed: other.use_n_observed,
        }
    }

    /// Move constructor analogue — Rust takes ownership by value, mirroring
    /// the C++ rvalue-reference constructor.
    pub fn move_construct(other: HyperLogLogPlusMinus) -> Self {
        HyperLogLogPlusMinus {
            p: other.p,
            m: other.m,
            registers: other.registers,
            n_observed: other.n_observed,
            sparse: other.sparse,
            sparse_list: other.sparse_list,
            use_n_observed: other.use_n_observed,
        }
    }

    /// Move-assignment analogue (`operator=(HyperLogLogPlusMinus&&)`).
    pub fn assign_from_move(&mut self, other: HyperLogLogPlusMinus) -> &mut Self {
        self.p = other.p;
        self.m = other.m;
        self.registers = other.registers;
        self.n_observed = other.n_observed;
        self.sparse = other.sparse;
        self.sparse_list = other.sparse_list;
        self.use_n_observed = other.use_n_observed;
        self
    }

    /// Copy-assignment analogue (`operator=(const HyperLogLogPlusMinus&)`).
    pub fn assign_from_copy(&mut self, other: &HyperLogLogPlusMinus) -> &mut Self {
        self.p = other.p;
        self.m = other.m;
        self.registers = other.registers.clone();
        self.n_observed = other.n_observed;
        self.sparse = other.sparse;
        self.sparse_list = other.sparse_list.clone();
        self.use_n_observed = other.use_n_observed;
        self
    }

    /// Observe a single 64-bit item. The item is mixed via
    /// `murmurhash3_finalizer` and either added to the sparse list (sparse
    /// mode) or used to update the corresponding register's rank (dense
    /// mode). Triggers a transition to dense mode when the sparse list
    /// would exceed `m/4` entries.
    pub fn insert(&mut self, item: u64) {
        self.n_observed += 1;
        let hash_value = murmurhash3_finalizer(item);

        if self.sparse && self.sparse_list.len() + 1 > self.m / 4 {
            self.switch_to_normal();
        }

        if self.sparse {
            let encoded = encode_hash_in_32bit(hash_value, P_PRIME, self.p);
            add_hash_to_sparse_set(&mut self.sparse_list, encoded, P_PRIME);
        } else {
            let idx = get_index_u64(hash_value, self.p) as usize;
            let rank = get_rank_u64(hash_value, self.p);
            if rank > self.registers[idx] {
                self.registers[idx] = rank;
            }
        }
    }

    /// Observe a batch of items by calling `insert` on each.
    pub fn insert_many(&mut self, items: &[u64]) {
        for &item in items {
            self.insert(item);
        }
    }

    /// Transition from sparse to dense representation, materialising the
    /// register array from the accumulated sparse list and clearing the
    /// sparse list. No-op if already in dense mode.
    fn switch_to_normal(&mut self) {
        if !self.sparse {
            return;
        }
        self.sparse = false;
        self.registers = vec![0u8; self.m];
        self.add_to_registers_from_sparse(&self.sparse_list.clone());
        self.sparse_list.clear();
    }

    /// Fold each entry in the given sparse list into our dense register
    /// array, decoding the index from the top `p` bits and taking the max
    /// of the current and encoded rank.
    fn add_to_registers_from_sparse(&mut self, sparse_list: &BTreeSet<u32>) {
        for &encoded in sparse_list {
            // Index extraction uses the upper p bits of the 32-bit encoded value
            let idx = (encoded >> (32 - self.p)) as usize;
            let rank_val = get_encoded_rank(encoded, P_PRIME, self.p);
            if rank_val > self.registers[idx] {
                self.registers[idx] = rank_val;
            }
        }
    }

    /// Merge another HLL into this one.
    pub fn merge(&mut self, other: &HyperLogLogPlusMinus) {
        assert_eq!(self.p, other.p, "precisions must be equal");
        if other.n_observed == 0 {
            return;
        }
        if self.n_observed == 0 {
            self.n_observed = other.n_observed;
            self.sparse = other.sparse;
            self.sparse_list = other.sparse_list.clone();
            self.registers = other.registers.clone();
            return;
        }

        self.n_observed += other.n_observed;
        if self.sparse && other.sparse {
            self.sparse_list.extend(other.sparse_list.iter());
        } else if other.sparse {
            self.add_to_registers_from_sparse(&other.sparse_list);
        } else if self.sparse {
            self.sparse = false;
            self.registers = other.registers.clone();
            self.add_to_registers_from_sparse(&self.sparse_list.clone());
            self.sparse_list.clear();
        } else {
            for i in 0..other.registers.len() {
                if other.registers[i] > self.registers[i] {
                    self.registers[i] = other.registers[i];
                }
            }
        }
    }

    /// Merge another HLL into this one, consuming it. Equivalent in effect
    /// to `merge`, but moves storage out of `other` to avoid extra copies
    /// (Rust analogue of the C++ `merge(HyperLogLogPlusMinus&&)`).
    pub fn merge_owned(&mut self, other: HyperLogLogPlusMinus) {
        assert_eq!(self.p, other.p, "precisions must be equal");
        if other.n_observed == 0 {
            return;
        }
        if self.n_observed == 0 {
            self.n_observed = other.n_observed;
            self.sparse = other.sparse;
            self.sparse_list = other.sparse_list;
            self.registers = other.registers;
            return;
        }

        self.n_observed += other.n_observed;
        if self.sparse && other.sparse {
            self.sparse_list.extend(other.sparse_list);
        } else if other.sparse {
            self.add_to_registers_from_sparse(&other.sparse_list);
        } else if self.sparse {
            self.sparse = false;
            self.registers = other.registers;
            self.add_to_registers_from_sparse(&self.sparse_list.clone());
            self.sparse_list.clear();
        } else {
            for i in 0..other.registers.len() {
                if other.registers[i] > self.registers[i] {
                    self.registers[i] = other.registers[i];
                }
            }
        }
    }

    /// Original Flajolet et al. HLL estimator with linear-counting fallback
    /// for small cardinalities (raw estimate ≤ 2.5m). If `use_sparse_precision`
    /// is true and we are in sparse mode, returns the linear-counting
    /// estimate over `m' = 2^p'` (the sparse-mode register count); otherwise
    /// promotes the sparse list to a temporary register array first.
    pub fn flajolet_cardinality(&self, use_sparse_precision: bool) -> u64 {
        let mut registers;
        let registers_ref = if self.sparse {
            if use_sparse_precision {
                return linear_counting(M_PRIME, M_PRIME - self.sparse_list.len() as u32).round()
                    as u64;
            }

            registers = vec![0u8; self.m];
            for &val in &self.sparse_list {
                let idx = get_index_u32(val, self.p) as usize;
                let rank_val = get_encoded_rank(val, P_PRIME, self.p);
                if rank_val > registers[idx] {
                    registers[idx] = rank_val;
                }
            }
            &registers
        } else {
            &self.registers
        };

        let mut est = calculate_raw_estimate(registers_ref);
        if est <= 2.5 * self.m as f64 {
            let v = count_zeros(registers_ref);
            if v > 0 {
                est = linear_counting(self.m as u32, v);
            }
        }

        if self.use_n_observed && (self.n_observed as f64) < est {
            self.n_observed
        } else {
            est.round() as u64
        }
    }

    /// HLL++ estimator of Heule et al. (2013). Uses empirical bias correction
    /// (`get_estimate_bias`) for raw estimates below `5m`, plus
    /// linear-counting fallback when the LC estimate is below the
    /// per-precision threshold. Only valid for `p <= 18`; falls back to
    /// the Ertl estimator otherwise.
    pub fn heule_cardinality(&self, correct_bias: bool) -> u64 {
        if self.p > 18 {
            return self.ertl_cardinality();
        }
        if self.sparse {
            return linear_counting(M_PRIME, M_PRIME - self.sparse_list.len() as u32).round()
                as u64;
        }

        let v = count_zeros(&self.registers);
        if v != 0 {
            let lc_estimate = linear_counting(self.m as u32, v).round() as u64;
            if lc_estimate <= THRESHOLD[(self.p - 4) as usize] as u64 {
                return lc_estimate;
            }
        }

        let mut est = calculate_raw_estimate(&self.registers);
        if correct_bias && est <= self.m as f64 * 5.0 {
            let bias = get_estimate_bias(est, self.p);
            assert!(est > bias);
            est -= bias;
        }

        if self.use_n_observed && (self.n_observed as f64) < est {
            self.n_observed
        } else {
            est.round() as u64
        }
    }

    /// Ertl cardinality estimator (default).
    pub fn ertl_cardinality(&self) -> u64 {
        let (q, m, c) = if self.sparse {
            let q = 64 - P_PRIME as usize;
            let m = M_PRIME as usize;
            let c = sparse_register_histogram(&self.sparse_list, P_PRIME, self.p, q);
            (q, m, c)
        } else {
            let q = 64 - self.p as usize;
            let m = self.m;
            let c = register_histogram(&self.registers, q);
            (q, m, c)
        };

        let mut est_denominator = m as f64 * tau(1.0 - c[q + 1] as f64 / m as f64);
        for k in (1..=q).rev() {
            est_denominator += c[k] as f64;
            est_denominator *= 0.5;
        }
        est_denominator += m as f64 * sigma(c[0] as f64 / m as f64);

        let m_sq_alpha_inf = (m as f64 / (2.0 * 2.0f64.ln())) * m as f64;
        let est = m_sq_alpha_inf / est_denominator;

        if self.use_n_observed && (self.n_observed as f64) < est {
            self.n_observed
        } else {
            est.round() as u64
        }
    }

    /// Default cardinality estimate.
    pub fn cardinality(&self) -> u64 {
        self.ertl_cardinality()
    }

    /// Alias for `cardinality()` — matches C++ `size()`.
    pub fn size(&self) -> u64 {
        self.cardinality()
    }

    /// Total number of items observed via `insert` (the multiset size, as
    /// opposed to the distinct-count estimate). Used as an upper bound on
    /// the cardinality when `use_n_observed` is set.
    pub fn n_observed(&self) -> u64 {
        self.n_observed
    }

    /// Reset the sketch to its empty state and return to sparse
    /// representation.
    pub fn reset(&mut self) {
        self.sparse = true;
        self.sparse_list.clear();
        self.registers.clear();
        self.n_observed = 0;
    }

    /// In-place merge by reference — analogue of C++
    /// `operator+=(const HyperLogLogPlusMinus&)`.
    pub fn add_assign_ref(&mut self, other: &HyperLogLogPlusMinus) -> &mut Self {
        self.merge(other);
        self
    }
}

impl Default for HyperLogLogPlusMinus {
    /// Default precision is 12 — matches the C++ default constructor.
    fn default() -> Self {
        Self::new(12)
    }
}

impl std::ops::AddAssign for HyperLogLogPlusMinus {
    /// `+=` merges another sketch into this one by consuming it.
    fn add_assign(&mut self, other: Self) {
        self.merge_owned(other);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hll_basic() {
        let mut hll = HyperLogLogPlusMinus::new(12);
        // Insert some values
        for i in 0..1000u64 {
            hll.insert(i);
        }
        let est = hll.cardinality();
        // Should be roughly 1000 (within ~5% for p=12)
        assert!(
            (900..=1100).contains(&est),
            "HLL estimate {} not close to 1000",
            est
        );
    }

    #[test]
    fn test_hll_merge() {
        let mut hll1 = HyperLogLogPlusMinus::new(12);
        let mut hll2 = HyperLogLogPlusMinus::new(12);

        for i in 0..500u64 {
            hll1.insert(i);
        }
        for i in 500..1000u64 {
            hll2.insert(i);
        }

        hll1.merge(&hll2);
        let est = hll1.cardinality();
        assert!(
            (900..=1100).contains(&est),
            "Merged HLL estimate {} not close to 1000",
            est
        );
    }

    #[test]
    fn test_hll_n_observed() {
        let mut hll = HyperLogLogPlusMinus::new(12);
        for i in 0..100u64 {
            hll.insert(i);
        }
        assert_eq!(hll.n_observed(), 100);
    }
}
