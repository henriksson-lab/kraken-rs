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

#[inline]
fn clz_u32(x: u32, max: u8) -> u8 {
    if x == 0 {
        max
    } else {
        x.leading_zeros() as u8
    }
}

#[inline]
fn clz_u64(x: u64, max: u8) -> u8 {
    if x == 0 {
        max
    } else {
        x.leading_zeros() as u8
    }
}

#[allow(dead_code)]
#[inline]
fn clz_p(x: u32, p: u8) -> u8 {
    (((x << 1) | 1) << (p - 1)).leading_zeros() as u8
}

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

#[allow(dead_code)]
#[inline]
fn extract_high_bits_u64(bits: u64, hi: u8) -> u64 {
    bits >> (64 - hi)
}

#[allow(dead_code)]
#[inline]
fn extract_high_bits_u32(bits: u32, hi: u8) -> u32 {
    bits >> (32 - hi)
}

#[inline]
fn get_index_u64(hash_value: u64, p: u8) -> u32 {
    (hash_value >> (64 - p)) as u32
}

#[allow(dead_code)]
#[inline]
fn get_index_u32(hash_value: u32, p: u8) -> u32 {
    hash_value >> (32 - p)
}

#[allow(dead_code)]
#[inline]
fn trailing_ones(p: u8) -> u64 {
    (1u64 << p) - 1
}

#[inline]
fn get_rank_u64(hash_value: u64, p: u8) -> u8 {
    let rank_bits = hash_value << p;
    clz_u64(rank_bits, 64 - p) + 1
}

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

#[inline]
fn get_rank_u32(hash_value: u32, p: u8) -> u8 {
    let rank_bits = hash_value << p;
    clz_u32(rank_bits, 32 - p) + 1
}

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

#[inline]
fn add_hash_to_sparse_set(set: &mut BTreeSet<u32>, val: u32, _p_prime: u8) {
    set.insert(val);
}

#[allow(dead_code)]
fn alpha(m: u32) -> f64 {
    match m {
        16 => 0.673,
        32 => 0.697,
        64 => 0.709,
        _ => 0.7213 / (1.0 + 1.079 / m as f64),
    }
}

#[allow(dead_code)]
fn linear_counting(m: u32, v: u32) -> f64 {
    assert!(v <= m, "number of v should not be greater than m");
    (m as f64) * ((m as f64) / (v as f64)).ln()
}

#[allow(dead_code)]
fn calculate_raw_estimate(registers: &[u8]) -> f64 {
    let m = registers.len();
    let mut inverse_sum: f64 = 0.0;
    for &r in registers {
        inverse_sum += 1.0 / (1u64 << r) as f64;
    }
    alpha(m as u32) * (m * m) as f64 / inverse_sum
}

#[allow(dead_code)]
fn count_zeros(registers: &[u8]) -> u32 {
    registers.iter().filter(|&&v| v == 0).count() as u32
}

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

fn register_histogram(registers: &[u8], q: usize) -> Vec<i32> {
    let mut c = vec![0i32; q + 2];
    for &r in registers {
        c[r as usize] += 1;
    }
    c
}

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
    pub fn new(precision: u8) -> Self {
        Self::with_sparse(precision, true)
    }

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

    pub fn insert_many(&mut self, items: &[u64]) {
        for &item in items {
            self.insert(item);
        }
    }

    fn switch_to_normal(&mut self) {
        if !self.sparse {
            return;
        }
        self.sparse = false;
        self.registers = vec![0u8; self.m];
        self.add_to_registers_from_sparse(&self.sparse_list.clone());
        self.sparse_list.clear();
    }

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

    pub fn size(&self) -> u64 {
        self.cardinality()
    }

    pub fn n_observed(&self) -> u64 {
        self.n_observed
    }

    pub fn reset(&mut self) {
        self.sparse = true;
        self.sparse_list.clear();
        self.registers.clear();
        self.n_observed = 0;
    }

    pub fn add_assign_ref(&mut self, other: &HyperLogLogPlusMinus) -> &mut Self {
        self.merge(other);
        self
    }
}

impl Default for HyperLogLogPlusMinus {
    fn default() -> Self {
        Self::new(12)
    }
}

impl std::ops::AddAssign for HyperLogLogPlusMinus {
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
