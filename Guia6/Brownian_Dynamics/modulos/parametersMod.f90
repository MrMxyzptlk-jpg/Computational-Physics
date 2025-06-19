MODULE parametersMod
    use precisionMod
    use constantsMod
    use functionsMod
    use mzranMod
    implicit none

 !   private     radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor
  !  private     positions_buffer, msd_counts, msd_count

    real(pr)                :: radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor
    real(pr), allocatable   :: positions_buffer(:,:,:)
    integer, allocatable    :: msd_counts(:)
    integer                 :: msd_count

    integer(int_medium)     :: cell_dim(3)
    character (len=15)      :: integrator
    character (len=6)       :: structure
    integer(int_large)      :: real_steps, measuring_jump, measuring_steps
    integer(int_large)      :: num_atoms, pair_corr_bins, max_correlation, MC_adjust_step
    real(pr)                :: conversion_factors(6), periodicity(3), lattice_constant, sigma, epsilon, dt, dtdt, Berendsen_time
    real(pr)                :: radius_cutoff, pair_corr_cutoff, dr, initial_Temp_Adim, density, mass, MC_delta
    logical                 :: transitory, save_transitory, do_pair_correlation, save_observables, measure

CONTAINS

subroutine initialize_parameters()

    if (structure == "random") then
        cell_dim = 1_int_medium
        density = num_atoms / product(cell_dim*lattice_constant)
    end if

    ! Define conversion factors to adimensionalize the variables
    conversion_factors(1) = sigma                           ! Distance      [Bohr = a₀]
    conversion_factors(2) = sigma*sqrt(mass/epsilon)        ! Time          [ℏ/Eh = tₐ]
    conversion_factors(3) = 1._pr                           ! Temperature   [epsilon / kB]
    conversion_factors(4) = epsilon                         ! Energy        [hartree = Eh]
    conversion_factors(5) = sqrt(epsilon/mass)              ! Velocity      [a₀ / tₐ]
    conversion_factors(6) = epsilon/(sigma*sigma*sigma)     ! Pressure      [hartree / bohr³]

    lattice_constant = lattice_constant/conversion_factors(1)
    initial_Temp_Adim = initial_Temp_Adim*conversion_factors(3)    ! Already non-dimensional!!
    periodicity = cell_dim*lattice_constant
    radius_cutoff = radius_cutoff/conversion_factors(1)

    dt = dt/conversion_factors(2)
    dtdt = dt*dt

    radius_cutoff_squared = radius_cutoff*radius_cutoff
    potential_cutoff = Lennard_Jones_potential(radius_cutoff_squared)

    transitory = .True.    ! Flag to avoid calculations and saving variables during the transitory steps

    if(do_pair_correlation) then
        pair_corr_cutoff_sqr = pair_corr_cutoff*pair_corr_cutoff
        dr = pair_corr_cutoff/real(pair_corr_bins,pr)
    end if

    call mzran_init() ! Initialize random number generator

    measure = .False.
    measuring_steps = real_steps/measuring_jump

end subroutine initialize_parameters

subroutine initialize_rest()

    Temp_factor = 2.0_pr / (3.0_pr * real(num_atoms,pr))
    Pressure_factor = 1._pr / (3._pr * product(cell_dim)*lattice_constant*lattice_constant*lattice_constant)

end subroutine initialize_rest



END MODULE parametersMod