MODULE parametersMod
    use precisionMod
    use constantsMod
    implicit none

    ! Internal variables
    real(pr)                :: volume, radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor
    real(pr)                :: dr, dtdt

    ! Observables' variables
    real(pr)                :: pair_corr_cutoff
    integer(int_large)      :: pair_corr_bins, max_correlation

    ! Physical variables
    real(pr)                :: conversion_factors(6), periodicity(3), lattice_constant, radius_cutoff, ref_Temp, density, mass, dt
    integer(int_large)      :: num_atoms
    integer                 :: cell_dim(3)
    character (len=6)       :: structure
    character (len=15)      :: integrator

    ! Main loops ranges
    integer(int_large)      :: real_steps, measuring_jump, measuring_steps

    ! Thermostats variables
    real(pr)                :: Berendsen_time

    ! Potential variables and internal variables
    real(pr)                :: sigma, epsilon               ! User variables
    real(pr)                :: sigma_sqr, Ewald_realFactor, twoPi_over_volume, eightPi_over_volume, halfSigma_sqr ! Internal variables (for Coulomb_Ewald)
    real(pr)                :: k_periodicity(3), Ewald_selfTerm
    integer                 :: num_kvec, kgrid(3), k_sqr_max    ! Internal variables (for Coulomb_Ewald)
    real(pr), allocatable   :: kfac(:)                          ! Internal variables (for Coulomb_Ewald)
    type :: kvector_data
        real(pr)    :: k_squared, kvec(3)
        real(pr)    :: k_factor
    end type kvector_data
    type(kvector_data), allocatable :: kvectors(:)   ! Internal variables (for Coulomb_Ewald)

    ! Monte Carlo specific variables
    integer(int_large)      :: MC_adjust_step
    real(pr)                :: MC_delta

    ! Brownian Dynamics specific variables and internal variables
    real(pr)                :: reduced_viscosity, viscosity                             ! User variables
    real(pr)                :: diffusion_coeff, reduced_viscosity_inv, brownian_stddev  ! Internal variables

    ! Tasks and measuring helper variables
    logical                 :: save_transitory, do_pair_correlation, save_observables, transitory, measure

    ! Potentials' pointers
    procedure(pot), pointer      :: potential => null()
    procedure(pot_func), pointer :: potential_function => null()

    abstract interface ! Intended to allow for the implementation of a different potential later on
        subroutine pot(particle_distance_squared, particle_separation, force_contribution, E_potential, pressure_virial &
            , potential_cutoff)
            use precisionMod
            real(pr), intent(in)               :: particle_distance_squared, particle_separation(3)
            real(pr), intent(out)              :: force_contribution(3)
            real(pr), intent(inout)            :: E_potential, pressure_virial
            real(pr), intent(in)               :: potential_cutoff
        end subroutine pot
        function pot_func(particle_distance_squared)
            use precisionMod
            real(pr), intent(in)    :: particle_distance_squared
            real(pr)                :: pot_func
        end function pot_func
    end interface

END MODULE parametersMod