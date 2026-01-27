! Created by  on 17/11/2023.

module ml_constants_mod
    use, intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, &
                                            real32, real64, real128
    implicit none
    private

    public :: r_ml, r_um, ml_nlev, ml_output_len, ml_input_len, &
              reference_pressure, poisson_constant

    integer, parameter :: r_ml = real32
    integer, parameter :: r_um = real64
    integer, parameter :: ml_nlev = 70
    integer, parameter :: ml_output_len = 70
    integer, parameter :: ml_input_len = 213
    real(kind=r_um), parameter :: reference_pressure = 100000
    real(kind=r_um), parameter :: poisson_constant = 0.2854
end module ml_constants_mod