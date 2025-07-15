! Module containing all global variables which are parsed or derived from parsed variables

MODULE variablesMod
    use precisionMod
    use constantsMod
    implicit none

!######### PHYSICAL variables ################
    ! Parsed
    character (len=6)       :: ensemble, structure
    real(pr)                :: lattice_constant, ref_Temp, density, mass
    integer                 :: cell_dim(3)                  ! Repetition of unit cell in each direction
    integer(int_large)      :: num_atoms                    ! Number of atoms
    real(pr)                :: reduced_viscosity, viscosity ! Only for Brownian Dynamics

    ! Internal
    real(pr)                :: conversion_factors(6), periodicity(3)
    real(pr)                :: reduced_viscosity_inv        ! Only for Brownian Dynamics

!######### CALCULATION variables #############
    ! Parsed
    character (len=11)      :: state            ! specifies if the starting configurations is read from file (fromFile) or started from scratch (fromScratch)
    real(pr)                :: radius_cutoff    ! radius cut-off for short range interactions
    real(pr)                :: dt               ! time step
    integer(int_large)      :: real_steps, transitory_steps, measuring_jump, measuring_steps ! Main loops ranges
    character (len=12)      :: summation        ! How pair interactions are computed. Short range: all-vs-all or linked-lists. Long range: Ewald
    character (len=12)      :: initial_velocities   ! Starting velocity distribution
    integer(int_large)      :: dim_linkCell(3)  ! Number of cells in each direction to sum over in linked-lists summation
    integer(int_large)      :: thermostat_steps

    ! Internal
    logical                 :: do_linkCell = .False.
    integer(int_large)      :: transitory_minIndex

!######### TASKS variables ###################
    ! Parsed
    logical                 :: save_positions, save_state, save_transitory, save_observables
    logical                 :: do_structure_factor, do_mean_sqr_displacement, do_pair_correlation
    integer(int_large)      :: Miller_index(3)

    ! Internal
    logical                 :: transitory, measure
    real(pr)                :: reciprocal_vec(3)

!######### APPROXIMATION variables ###########
    ! Parsed
    character (len=15)      :: integrator, interactions
    real(pr)                :: sigma, delta   ! Variables for the potential defining the interactions
    integer                 :: kgrid(3)         ! Partition of reciprocal space ('Ewald' summation)
    real(pr)                :: MC_acceptance_rate   ! ('Monte-Carlo' integrator)
    real(pr)                :: MC_delta         ! ('Monte-Carlo' integrator)

    ! Internal (Ewald summation)
    real(pr)                :: sigma_sqr, Ewald_realFactor, twoPi_over_volume, eightPi_over_volume, halfSigma_sqr, Ewald_selfTerm ! To avoid recalculation
    real(pr)                :: k_periodicity(3) ! reciprocal lattice periodicity
    integer                 :: num_kvec         ! Number of allowed reciprocal vectors
    integer                 :: k_sqr_max        ! Maximum reciprocal vector length allowed
    real(pr), allocatable   :: kfac(:)          ! Factors appearing in 'Ewald' summation (avoid recalculation)
    complex(pr), allocatable    :: reciprocal_charges(:)    ! Fourier transform of the periodic charge distribution (avoid recalculation in MC)
   ! type :: kvector_data       ! UNUSED
   !     real(pr)    :: k_squared, kvec(3),k_factor
   ! end type kvector_data
   ! type(kvector_data), allocatable :: kvectors(:)

!######### THERMOSTAT variables ##############
    ! Parsed
    character (len=12)      :: thermostat_type
    real(pr)                :: Berendsen_time

    ! Internal
    integer(int_large)      :: MC_accepted

!######### MSD variables #####################
    ! Parsed
    integer(int_large)      :: max_correlation

!######### PAIR_CORRELATION variables ########
    ! Parsed
    real(pr)                :: pair_corr_cutoff     ! Cut-off radius (must be smaller than the minimum periodicity length)
    integer(int_large)      :: pair_corr_bins       ! Number of bins for the histogram

!######### OTHER INTERNAL VARIABLES ##########

    ! Common variables
    real(pr)                :: volume, Temp_factor, Pressure_factor
    real(pr)                :: dr, dtdt
    real(pr)                :: radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff


    ! Brownian Dynamics specific variables and internal variables
    real(pr)                :: diffusion_coeff, brownian_stddev  ! Internal variables

END MODULE variablesMod