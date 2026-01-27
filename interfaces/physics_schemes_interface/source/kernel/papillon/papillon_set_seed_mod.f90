module papillon_set_seed_mod
    use papillon_config_mod, only: random_seed_size, random_seed_value
    implicit none
contains
    subroutine random_seed_set_fixed()
        call random_seed(size=random_seed_size)
        call random_seed(put=random_seed_value)
    end subroutine random_seed_set_fixed
end module papillon_set_seed_mod