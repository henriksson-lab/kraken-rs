use std::collections::BTreeSet;

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

#[inline]
fn clz_u32(x: u32, max: u8) -> u8 {
    if x == 0 { max } else { x.leading_zeros() as u8 }
}

#[inline]
fn clz_u64(x: u64, max: u8) -> u8 {
    if x == 0 { max } else { x.leading_zeros() as u8 }
}

#[inline]
fn get_index_u64(hash_value: u64, p: u8) -> u32 {
    (hash_value >> (64 - p)) as u32
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
    let idx = ((hash_value >> (64 - p_prime)) as u32) << (32 - p_prime);
    if (idx << p) == 0 {
        let additional_rank = get_rank_u64(hash_value, p_prime);
        idx | ((additional_rank as u32) << 1) | 1
    } else {
        idx
    }
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

fn register_histogram(registers: &[u8], q: usize) -> Vec<i32> {
    let mut c = vec![0i32; q + 2];
    for &r in registers {
        c[r as usize] += 1;
    }
    c
}

fn sparse_register_histogram(sparse_list: &BTreeSet<u32>, p_prime: u8, p: u8, q: usize) -> Vec<i32> {
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

/// Ertl's tau correction for high-value registers.
fn tau(mut x: f64) -> f64 {
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
        assert!((4..=18).contains(&precision),
            "precision must be between 4 and 18");
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

    pub fn insert(&mut self, item: u64) {
        self.n_observed += 1;
        let hash_value = murmurhash3_finalizer(item);

        if self.sparse && self.sparse_list.len() + 1 > self.m / 4 {
            self.switch_to_normal();
        }

        if self.sparse {
            let encoded = encode_hash_in_32bit(hash_value, P_PRIME, self.p);
            self.sparse_list.insert(encoded);
        } else {
            let idx = get_index_u64(hash_value, self.p) as usize;
            let rank = get_rank_u64(hash_value, self.p);
            if rank > self.registers[idx] {
                self.registers[idx] = rank;
            }
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
}

impl Default for HyperLogLogPlusMinus {
    fn default() -> Self {
        Self::new(12)
    }
}

impl std::ops::AddAssign for HyperLogLogPlusMinus {
    fn add_assign(&mut self, other: Self) {
        self.merge(&other);
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
        assert!(est >= 900 && est <= 1100,
            "HLL estimate {} not close to 1000", est);
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
        assert!(est >= 900 && est <= 1100,
            "Merged HLL estimate {} not close to 1000", est);
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
