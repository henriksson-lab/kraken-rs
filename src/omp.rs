pub type OmpLockT = i32;

pub fn omp_get_thread_num() -> i32 {
    0
}

pub fn omp_get_max_threads() -> i32 {
    1
}

pub fn omp_set_num_threads(_num: i32) {}

pub fn omp_init_lock(_lock: &mut OmpLockT) {}

pub fn omp_destroy_lock(_lock: &mut OmpLockT) {}

pub fn omp_set_lock(_lock: &mut OmpLockT) {}

pub fn omp_unset_lock(_lock: &mut OmpLockT) {}

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
