!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Kernel for the PAPILLON stochastic physics scheme.
!>
module papillon_kernel_mod

  use argument_mod,           only : arg_type,          &
                                     GH_FIELD, GH_REAL, &
                                     GH_SCALAR, &
                                     GH_READ, GH_WRITE, &
                                     CELL_COLUMN,       &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2, &
                                     ANY_DISCONTINUOUS_SPACE_4
  use constants_mod,          only : i_def, i_um, r_def, r_um
  use fs_continuity_mod,      only : W3, Wtheta
  use kernel_mod,             only : kernel_type
  use papillon_config_mod,    only: lat_scale_factor, lon_scale_factor, &
                                 radius_scale_factor, height_scale_factor, time_scale_factor, &
                                 use_land_fraction, use_orog, use_sd_orog, noise_scale_factor
  use extrusion_config_mod,      only: planet_radius
  use ml_constants_mod,       only: r_ml, reference_pressure, poisson_constant, &
          ml_nlev, ml_input_len
  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: papillon_kernel_type
    private
    type(arg_type) :: meta_args(15) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dtheta_papillon
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! noise
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! theta_star
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sd_orog
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! time_seconds
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),  & ! latitude
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2),  & ! longitude
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! noise_centre_x
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! noise_centre_y
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            & ! noise_centre_z
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                             & ! noise_centre_t
            /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: papillon_code
  end type papillon_kernel_type

  public :: papillon_code

contains

  !> @brief Implementation of the PAPILLON stochastic physics scheme
  !> @details The PAPILLON stochastic physics scheme is a
  !>           stochastic physics parametrization that diagnoses
  !>           subgrid variability with a machine learning model
  !>           and perturbs the inputs to the other parameterizations
  !> @param[in]     nlayers          Number of layers
  !> @param[out]    dtheta_papillon  Papillon temperature increment
  !> @param[in,out] noise            Papillon noise field
  !> @param[in]     theta_star       Potential temperature predictor after advection
  !> @param[in]     m_v              Vapour mixing ratio after advection
  !> @param[in]     exner_in_wth     Exner pressure field in wth space
  !> @param[in]     height_wth       Altitude on wth levels
  !> @param[in]     sd_orog          Standard deviation of subgrid orography
  !> @param[in]     tile_fraction    Fraction of surface that is land
  !> @param[in]     time_seconds     Simulation time, used so noise varies in space and time
  !> @param[in]     latitude         Latitude of the cell
  !> @param[in]     longitude        Longitude of the cell
  !> @param[in]     noise_centre_x   Centre of the noise field in the x dimension
  !> @param[in]     noise_centre_y   Centre of the noise field in the y dimension
  !> @param[in]     noise_centre_z   Centre of the noise field in the z dimension
  !> @param[in]     noise_centre_t   Centre of the noise field in the t dimension
  !> @param[in]     ndf_wth          Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth         Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth          Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_2d           Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d          Number of unique DOFs for 2D fields
  !> @param[in]     map_2d           Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile         Number of DOFs per cell for tiles
  !> @param[in]     undf_tile        Number of total DOFs for tiles
  !> @param[in]     map_tile         Dofmap for cell for surface tiles
  subroutine papillon_code(nlayers,         &
                          dtheta_papillon, &
                          noise,           &
                          theta_star,      &
                          m_v,             &
                          exner_in_wth,    &
                          height_wth,      &
                          sd_orog,         &
                          tile_fraction,   &
                          time_seconds,    &
                          latitude,        &
                          longitude,       &
                          noise_centre_x,  &
                          noise_centre_y,  &
                          noise_centre_z,  &
                          noise_centre_t,  &
                          ndf_wth,         &
                          undf_wth,        &
                          map_wth,         &
                          ndf_2d,          &
                          undf_2d,         &
                          map_2d,          &
                          ndf_tile,        &
                          undf_tile,       &
                          map_tile         &
          )
    !---------------------------------------
    ! Lfric modules
    !---------------------------------------
    use plausible_alternative_mod, only: get_alternative_temperature
    use planet_constants_mod, only: p_zero, kappa
    use simplex_noise_mod, only: snoise4d
    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use nlsizes_namelist_mod, only: row_length, rows
    ! TODO: Remove dependence on UM modules
    
    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, ndf_2d
    integer(kind=i_def), intent(in) :: undf_wth,undf_2d, ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile)

    integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
    integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d

    real(kind=r_def), dimension(undf_wth), intent(inout)  :: dtheta_papillon, noise

    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_star,   &
                                                           m_v,          &
                                                           exner_in_wth, &
                                                           height_wth
    
    real(kind=r_def), dimension(undf_2d),  intent(in)  :: sd_orog, latitude, longitude
    real(kind=r_def), dimension(undf_tile),  intent(in)  :: tile_fraction
    real(kind=r_def), intent(in) :: time_seconds
    real(r_def), intent(in) :: noise_centre_x, noise_centre_y, noise_centre_z, noise_centre_t

    ! Local variables for the kernel
    integer(kind=i_def) :: i, k

    real(r_def), dimension(nlayers) :: theta_papillon, q_papillon, &
    theta_inc, q_inc, p_theta_levels, &
    papillon_theta_inc, papillon_q_inc, &
    temperature_papillon, alternative_theta_papillon
    real(r_um), dimension(nlayers) :: noise_papillon
    real(r_def), dimension(4) :: noise_loc
!    real(r_def), dimension(row_length,rows) ::   p_star

    ! ML local variables
!    integer :: i
    real(kind=r_ml) :: inputs(ml_input_len)
    real(kind=r_def) :: alternative_temperature(ml_nlev)
    real(kind=r_def) :: sin_latitude, cos_latitude, sin_longitude, cos_longitude, &
                        noise_radius

    !-----------------------------------------------------------------------
    ! Initialise variables required from input fields
    !-----------------------------------------------------------------------
    do k = 1, nlayers
      theta_papillon(k) = theta_star(map_wth(1) + k)
      q_papillon(k) = m_v(map_wth(1) + k)
      p_theta_levels(k) = p_zero*(exner_in_wth(map_wth(1) + k))**(1.0_r_def/kappa)
    end do

    ! Initialise increments and output fields to zero
    theta_inc(:) = 0.0_r_def
    q_inc(:) = 0.0_r_def
    ! Generate the noise field
    noise_loc = [0.0_r_def, 0.0_r_def, 0.0_r_def, time_scale_factor*time_seconds - noise_centre_t]
    sin_latitude = sin(latitude(map_2d(1)))
    cos_latitude = cos(latitude(map_2d(1)))
    sin_longitude = sin(longitude(map_2d(1)))
    cos_longitude = cos(longitude(map_2d(1)))
    do i=1, nlayers
      noise_radius = (radius_scale_factor*planet_radius + height_scale_factor*height_wth(map_wth(1)+i))
      noise_loc(1) = noise_radius*sin_latitude*cos_longitude - noise_centre_x
      noise_loc(2) = noise_radius*sin_latitude*sin_longitude - noise_centre_y
      noise_loc(3) = noise_radius*cos_latitude - noise_centre_z
      noise_papillon(i) = noise_scale_factor*snoise4d(noise_loc)
    end do
    !-----------------------------------------------------------------------
    ! Call the ML model
    !-----------------------------------------------------------------------
!    inputs = [1.0_r_ml, 2500._r_ml, 300._r_ml, 0.64_r_ml, 0.9_r_ml, 0.91_r_ml]
!    inputs = [1.0_r_ml, 0.12_r_ml, 0.02_r_ml, 0.64_r_ml, 0.9_r_ml, 0.91_r_ml]
!    inputs = [0.0_r_ml, 0.0_r_ml, 0.0_r_ml, 0.17_r_ml, 0.92_r_ml, 0.75_r_ml]
    temperature_papillon = theta_papillon * (p_theta_levels / reference_pressure)**poisson_constant
    ! set ancillary inputs according to config
    if (use_land_fraction) then
      inputs(211) = real(tile_fraction(1), r_ml)
    else
      inputs(211) = 0_r_ml
    end if
    if (use_orog) then
      ! TODO: replace height of first model level with orography
      inputs(212) = real(height_wth(map_wth(1)+1), r_ml)
    else
      inputs(212) = 0_r_ml
    end if
    if (use_sd_orog) then
      inputs(213) = real(sd_orog(1), r_ml)
    else
      inputs(213) = 0_r_ml
    end if
    inputs(1:70) = real(q_papillon(:), r_ml)
    inputs(71:140) = real(p_theta_levels(:), r_ml)
    inputs(141:210) = real(temperature_papillon(:), r_ml)
!    print*, 'temperature size:', &
!            'theta conv(1,1,:)', size(theta_papillon(:)),&
!            'theta star', size(theta_star),&
!            'nlayers', nlayers,&
!            'temperature_papillon', size(temperature_papillon)
!        print*, 'inputs', inputs
    call get_alternative_temperature(inputs, noise_papillon, theta_papillon(:), alternative_temperature)
    alternative_theta_papillon(:) = &
            alternative_temperature
    ! Write out final increment
    do k = 1, nlayers
      dtheta_papillon(map_wth(1)+k) = alternative_theta_papillon(k) - theta_papillon(k)
      noise(map_wth(1)+k) = noise_papillon(k)
    end do
!    print*, 'plausible alternative'
!    do i=1, nlayers
!      print*, i, &
!              'mean T', real(temperature_papillon(i), r_ml), &
!              'pa T', real(alternative_temperature(i), r_ml), &
!              'mean theta', real(theta_papillon(i), r_ml), &
!              'pa theta', real(alternative_theta_papillon(i),r_ml)
!    end do
!        print*, 'tile_fraction', tile_fraction
!        print*, 'height_wth(1)', height_wth(1)
!        print*, 'sd_orog', sd_orog
!    print*, 'undf_2d', undf_2d
!    print*, 'size and shape of land area fraction', size(tile_fraction), shape(tile_fraction)


  end subroutine papillon_code

end module papillon_kernel_mod
