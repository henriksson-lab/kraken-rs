//! Single-threaded shims for the OpenMP runtime calls used by ported C++
//! code. Real parallelism in the Rust port is done with `rayon` /
//! `parking_lot`; these stubs exist so verbatim ports of OpenMP-aware
//! routines compile without restructuring.

/// Placeholder for `omp_lock_t`. The shim is single-threaded so the value
/// is never actually consulted.
pub type OmpLockT = i32;

/// Always returns 0 — there is only one (logical) thread in the shim.
pub fn omp_get_thread_num() -> i32 {
    0
}

/// Always returns 1 — the shim reports a single worker thread.
pub fn omp_get_max_threads() -> i32 {
    1
}

/// No-op: the shim ignores requests to change the OpenMP thread count.
pub fn omp_set_num_threads(_num: i32) {}

/// No-op: locks are unnecessary in single-threaded execution.
pub fn omp_init_lock(_lock: &mut OmpLockT) {}

/// No-op companion to [`omp_init_lock`].
pub fn omp_destroy_lock(_lock: &mut OmpLockT) {}

/// No-op: acquire is trivially successful in single-threaded mode.
pub fn omp_set_lock(_lock: &mut OmpLockT) {}

/// No-op: release pairs with [`omp_set_lock`].
pub fn omp_unset_lock(_lock: &mut OmpLockT) {}

/// Always returns 1 — the lock is reported as successfully acquired.
pub fn omp_test_lock(_lock: &mut OmpLockT) -> i32 {
    1
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_omp_hack_defaults() {
        let mut lock = 0;
        assert_eq!(omp_get_thread_num(), 0);
        assert_eq!(omp_get_max_threads(), 1);
        omp_set_num_threads(8);
        omp_init_lock(&mut lock);
        omp_set_lock(&mut lock);
        omp_unset_lock(&mut lock);
        omp_destroy_lock(&mut lock);
        assert_eq!(omp_test_lock(&mut lock), 1);
    }
}
