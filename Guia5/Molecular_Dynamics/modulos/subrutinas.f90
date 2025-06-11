MODULE subrutinas
    use precision
    use constantes
    use formats
    use funciones
    use mzranmod
    use omp_lib
    implicit none

    private     radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor
    private     positions_buffer, msd_counts, msd_count
    real(pr)                :: radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor
    real(pr), allocatable   :: positions_buffer(:,:,:)
    integer, allocatable    :: msd_counts(:)
    integer                 :: msd_count

    integer(int_medium)     :: cell_dim(3)
    character (len=6)       :: structure
    integer(int_large)      :: real_steps, measuring_jump, measuring_steps
    integer(int_large)      :: num_atoms, pair_corr_bins, max_correlation, MC_adjust_step
    real(pr)                :: conversion_factors(6), periodicity(3), lattice_constant, sigma, epsilon, dt, dtdt, Berendsen_time
    real(pr)                :: radius_cutoff, pair_corr_cutoff, dr, initial_Temp_Adim, density, mass, MC_delta
    logical                 :: transitory, save_transitory, do_pair_correlation, measure
    procedure(pot), pointer :: potential => null()
    procedure(pot_func), pointer :: potential_function => null()

    abstract interface ! Intended to allow for the implementation of a different potential later on
        subroutine pot(particle_distance_squared, particle_separation, force_contribution, E_potential, pressure_virial &
            , potential_cutoff)
            use precision
            real(pr), intent(in)               :: particle_distance_squared, particle_separation(3)
            real(pr), intent(out)              :: force_contribution(3)
            real(pr), intent(inout)            :: E_potential, pressure_virial
            real(pr), intent(in)               :: potential_cutoff
        end subroutine pot
        function pot_func(particle_distance_squared)
            use precision
            real(pr), intent(in)    :: particle_distance_squared
            real(pr)                :: pot_func
        end function pot_func
    end interface

contains

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

subroutine thermostat_rescale(velocities)
    real(pr), allocatable, intent(inout)    :: velocities(:,:)
    real(pr)                                :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities)/(3.0_pr * real(num_atoms,pr))
    scaling_factor = sqrt( initial_Temp_Adim / instant_Temp )

    velocities = velocities*scaling_factor

end subroutine thermostat_rescale

subroutine thermostat_Berendsen(velocities)
    real(pr), allocatable, intent(inout)    :: velocities(:,:)
    real(pr)                                :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities)/(3.0_pr * real(num_atoms,pr))

    scaling_factor = sqrt( 1._pr + dt/Berendsen_time*(initial_Temp_Adim / instant_Temp -1._pr))

    velocities = velocities*scaling_factor

end subroutine thermostat_Berendsen

subroutine get_forces_allVSall(positions, forces, E_potential, pressure_virial, pair_corr)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer(int_huge)       :: i, j

    forces = 0._pr
    if (measure) then; E_potential = 0.0; pressure_virial = 0.0 ; end if

    !$omp parallel private(j, i) &
    !$omp shared(positions, num_atoms) &
    !$omp reduction(+: forces, E_potential, pressure_virial, pair_corr)

        !$omp do schedule(dynamic)
        do i=1,num_atoms-1
            do j = i+1, num_atoms
                call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential, &
                    pressure_virial, pair_corr)
            end do
        end do
        !$omp end do

    !$omp end parallel

end subroutine get_forces_allVSall

subroutine get_force_contribution(particle1_position, particle2_position, particle1_forces, particle2_forces, E_potential &
    , pressure_virial, pair_corr)
    real(pr), dimension(3), intent(in)      :: particle1_position, particle2_position
    real(pr), dimension(3), intent(out)     :: particle1_forces, particle2_forces
    real(pr), intent(out)                   :: E_potential, pressure_virial
    real(pr), intent(inout)                 :: pair_corr(:)
    real(pr), dimension(3)                  :: particle_separation, force_contribution
    real(pr)                                :: particle_distance_squared

    particle_separation = particle1_position - particle2_position ! Separation vector
    particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
    particle_distance_squared = sum(particle_separation*particle_separation)

    if (particle_distance_squared <= radius_cutoff_squared) then
        call potential(particle_distance_squared, particle_separation, force_contribution, E_potential &
            , pressure_virial, potential_cutoff)
        particle1_forces = particle1_forces + force_contribution
        particle2_forces = particle2_forces - force_contribution
    endif

    if (do_pair_correlation .and. .not. transitory) call update_pair_correlation(particle_distance_squared, pair_corr)

end subroutine get_force_contribution

subroutine Lennard_Jones(particle_distance_squared, particle_separation,  force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_squared, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff
    real(pr)                   :: r2inv, r6inv, force_magnitude

    r2inv = 1._pr/particle_distance_squared
    r6inv = r2inv*r2inv*r2inv
    force_magnitude = 48._pr*r2inv*r6inv*(r6inv-0.5_pr)
    force_contribution = force_magnitude*particle_separation

    if (measure) then
        E_potential = E_potential + 4._pr*r6inv*(r6inv-1._pr) - potential_cutoff
        pressure_virial = pressure_virial + particle_distance_squared*force_magnitude
    end if

end subroutine Lennard_Jones

subroutine get_distance_squared(particle1, particle2, distance_squared)
    real(pr), intent(in)    :: particle1(3), particle2(3)
    real(pr)                :: distance_squared
    real(pr)                :: particle_separation(3)

    particle_separation = particle1 - particle2 ! Separation vector
    particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
    distance_squared = sum(particle_separation*particle_separation)

end subroutine get_distance_squared

subroutine update_positions_random(positions, E_potential, N_accepted)
    real(pr), dimension(:,:), intent(inout) :: positions
    real(pr), intent(inout)                 :: E_potential
    integer, intent(inout)                  :: N_accepted
    real(pr)                                :: random_displacement(3)
    real(pr)                                :: old_position(3), E_potential_old, E_potential_new, dE
    integer                                 :: random_particle_id, i, j

    do j = 1, num_atoms
        ! Pick a random particle
        random_particle_id = int(rmzran()*num_atoms) + 1
        old_position = positions(:,random_particle_id)

        ! Compute potential energy contribution
        call get_E_potential_contribution(positions, random_particle_id, E_potential_old)

        ! Propose a displacement
        random_displacement = (/(MC_delta*(rmzran() - 0.5d0), i = 1, 3)/)
        positions(:,random_particle_id) = old_position + random_displacement

        ! Apply periodic boundary conditions
        positions(:,random_particle_id) = modulo(positions(:,random_particle_id), periodicity(:))

        ! Compute new potential energy contribution
        call get_E_potential_contribution(positions, random_particle_id, E_potential_new)

        dE = E_potential_new - E_potential_old

        ! Metropolis criterion
        if (dE <= 0._pr ) then
            if (measure) E_potential = E_potential_new
            N_accepted = N_accepted + 1
            return
        else if ( rmzran() < exp(-dE / initial_Temp_Adim)) then
            if (measure) E_potential = E_potential_new
            N_accepted = N_accepted + 1
            return
        else
            ! Revert move if trial not accepted
            positions(:,random_particle_id) = old_position
        end if
    end do

end subroutine update_positions_random

subroutine get_E_potential_contribution(positions, random_particle_id, dE)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_squared
    integer                 :: random_particle_id, i

    dE = 0._pr
    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_squared)
            if (particle_distance_squared <= radius_cutoff_squared) then
                dE = dE + potential_function(particle_distance_squared) ! The term "- potential_cutoff" is irrelevant to the change in potential energy
            end if
        end if
    end do

end subroutine get_E_potential_contribution

subroutine update_random_step(N_accepted)
    integer, intent(inout)  :: N_accepted

    if (real(N_accepted,pr)/real(MC_adjust_step,pr) > 0.5_pr) then
        MC_delta = MC_delta*1.05_pr
    else
        MC_delta = MC_delta*0.95_pr
    endif
    N_accepted = 0

end subroutine update_random_step

subroutine update_positions_velVer(positions, velocities, forces)
    real(pr), dimension(:,:), intent(inout)    :: positions, velocities, forces

    positions = positions + velocities*dt + forces*dtdt*0.5_pr

    ! Apply periodic boundary conditions
    !positions = mod(positions, spread(periodicity, dim=2, ncopies=size(positions,2)))
    positions(1,:) = modulo(positions(1,:), periodicity(1))
    positions(2,:) = modulo(positions(2,:), periodicity(2))
    positions(3,:) = modulo(positions(3,:), periodicity(3))

end subroutine update_positions_velVer

subroutine update_velocities_velVer(velocities, forces, previous_forces)
    real(pr), dimension(:,:), intent(inout)    :: velocities, forces, previous_forces

    velocities = velocities + dt * 0.5_pr*(previous_forces + forces)

end subroutine update_velocities_velVer

subroutine get_pair_correlation_allVSall(positions, pair_corr)
    real(pr), intent(in)        :: positions(:,:)
    real(pr), intent(inout)     :: pair_corr(:)
    real(pr)                    :: particle_distance_squared
    integer(int_large)          :: i,j

    if (measure) then
        !$omp parallel private(j, i,  particle_distance_squared) &
        !$omp shared(positions, num_atoms) &
        !$omp reduction(+: pair_corr)

            !$omp do schedule(dynamic)
            do i=1,num_atoms-1
                do j = i+1, num_atoms
                    call get_distance_squared(positions(:,i), positions(:,j), particle_distance_squared)
                    call update_pair_correlation(particle_distance_squared, pair_corr)
                end do
            end do
            !$omp end do

        !$omp end parallel
    end if

end subroutine get_pair_correlation_allVSall

subroutine update_pair_correlation(particle_distance_squared, pair_corr)
    real(pr), intent(in)        :: particle_distance_squared
    real(pr), intent(inout)     :: pair_corr(:)
    integer(int_large)          :: radius_bin

    if (measure) then
        if (particle_distance_squared <= pair_corr_cutoff_sqr) then
            radius_bin = int(sqrt(particle_distance_squared) / dr) + 1
            pair_corr(radius_bin) = pair_corr(radius_bin) + 1._pr
        end if
    end if

end subroutine update_pair_correlation

subroutine normalize_pair_correlation(pair_corr)
    real(pr), intent(inout)     :: pair_corr(:)
    integer                     :: i
    real(pr)                    :: r_lower, r_upper, shell_vol, ideal_pair_corr

    do i = 1, pair_corr_bins
        r_lower = real(i-1,pr) * dr
        r_upper = r_lower + dr
        shell_vol = (4.0_pr*pi / 3.0_pr) * (r_upper**3 - r_lower**3)
        ideal_pair_corr = density * shell_vol * real(num_atoms,pr)
        pair_corr(i) = pair_corr(i) * 2._pr / ideal_pair_corr / real(measuring_steps,pr)
    end do

end subroutine normalize_pair_correlation

subroutine get_observables(velocities, E_kinetic, Pressure, Temperature)
    real(pr), intent(in)    :: velocities(:,:)
    real(pr), intent(inout) :: Pressure  ! Comes in as Pressure_virial
    real(pr), intent(out)   :: E_kinetic, Temperature

    E_kinetic = 0.5_pr*sum(velocities*velocities)
    Temperature = Temp_factor*E_kinetic
    Pressure = density*Temperature + Pressure*Pressure_factor

end subroutine get_observables

subroutine get_stats(measurements, variance, stddev, average)
    real(pr), intent(in)                :: measurements(:)
    real(pr), optional, intent(out)     :: variance, stddev, average
    real(pr)                            :: avg, sqr_avg

    if (present(average) .or. present(variance) .or. present(stddev)) then
        avg = sum(measurements)/size(measurements)
        sqr_avg = sum(measurements*measurements)/size(measurements)
    end if

    if  (present(average))   average  = avg
    if  (present(variance))  variance = sqr_avg - avg*avg
    if  (present(stddev))     stddev   = sqrt(sqr_avg - avg*avg)

end subroutine get_stats

subroutine get_reciprocal_vec(Miller_index, reciprocal_vec)
    integer, intent(in)     :: Miller_index(3)
    real(pr), intent(out)   :: reciprocal_vec(3)
    real(pr)                :: reciprocal_basis(3,3), reciprocal_latticeParam
    integer                 :: i

    select case (structure)
        case ("FCC")
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/-0.5_pr,  0.5_pr,  0.5_pr/)  ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/ 0.5_pr, -0.5_pr,  0.5_pr/)  ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/ 0.5_pr,  0.5_pr, -0.5_pr/)  ! reciprocal lattice basis vector 3
        case ("BCC")
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/0.0_pr, 0.5_pr, 0.5_pr/)      ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/0.5_pr, 0.0_pr, 0.5_pr/)      ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/0.5_pr, 0.5_pr, 0.0_pr/)      ! reciprocal lattice basis vector 3
        case ("random")
            print*, "Random structure selected -> Simple cubic basis set will be employed for structure factor calculation"
            reciprocal_latticeParam = 2._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/1._pr, 0._pr, 0._pr/)         ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/0._pr, 1._pr, 0._pr/)         ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/0._pr, 0._pr, 1._pr/)         ! reciprocal lattice basis vector 3
        case default
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/-0.5_pr,  0.5_pr,  0.5_pr/)   ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/ 0.5_pr, -0.5_pr,  0.5_pr/)   ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/ 0.5_pr,  0.5_pr, -0.5_pr/)   ! reciprocal lattice basis vector 3
    end select

    reciprocal_basis = reciprocal_basis * reciprocal_latticeParam

    reciprocal_vec = 0._pr
    do i =1, 3
        reciprocal_vec = reciprocal_vec + real(Miller_index(i),pr)*reciprocal_basis(:,i)
    end do

end subroutine get_reciprocal_vec

subroutine get_structure_factor(positions, structure_factor, reciprocal_vec)
    real(pr), intent(in)        :: positions(:,:), reciprocal_vec(3)
    real(pr), intent(inout)     :: structure_factor
    complex(pr)                 :: summation
    integer                     :: i

    summation = (0._pr,0._pr)

    do i = 1, num_atoms
        summation = summation + exp(CMPLX(0._pr,sum(reciprocal_vec*positions(:,i)), pr))
    end do

    structure_factor = abs(summation)*abs(summation)/(num_atoms*num_atoms)

end subroutine get_structure_factor

subroutine initialize_msd(msd)
    real(pr), allocatable, intent(out)  :: msd(:)

    allocate(positions_buffer(3,num_atoms, 0:max_correlation))
    allocate(msd_counts(0:max_correlation))
    allocate(msd(0:max_correlation))

    positions_buffer = 0._pr
    msd_counts = 0
    msd = 0._pr

end subroutine initialize_msd

subroutine update_msd(positions, msd)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(inout) :: msd(0:)
    integer                 :: j, ilast, inext
    real(pr)                :: displacements(3,num_atoms)

    msd_count = msd_count + 1
    ilast = mod(msd_count, max_correlation + 1)

    positions_buffer (:,:,ilast)= positions(:,:)

    do j = 0, min(max_correlation, msd_count-1)
        inext = mod(msd_count - j, max_correlation + 1)

        displacements = positions_buffer(:,:,ilast) - positions_buffer(:,:,inext)

        msd (j) = msd(j) + sum(displacements*displacements)
        msd_counts(j) = msd_counts(j) + 1
    end do

end subroutine update_msd

subroutine normalize_msd(msd)
    real(pr), intent(out)   :: msd(0:)
    integer                 :: j

    do j = 0, max_correlation
        if (msd_counts(j) > 0) then
            msd(j) = msd(j) / (real(msd_counts(j), pr) * real(num_atoms, pr))
        else
            msd(j) = 0._pr
        end if
    end do

end subroutine normalize_msd

subroutine check_measuring(index)
    integer, intent(in)         :: index

    if (mod(index,measuring_jump) == 0 ) then
        measure = .true.
    else
        measure = .false.
    end if

end subroutine check_measuring

subroutine gasdev_v(harvest)
    real(pr), dimension(:), intent(out) :: harvest
    real(pr), allocatable, save         :: g(:)
    logical, save                       :: gaus_stored = .false.
    integer, save                       :: last_allocated = 0
    real(pr)                            :: rsq, v1, v2
    integer                             :: i, n

    n = size(harvest)

    if (n /= last_allocated) then
        if (last_allocated > 0) deallocate(g)
        allocate(g(n))
        last_allocated = n
        gaus_stored = .false.
    end if

    if (gaus_stored) then
        harvest = g
        gaus_stored = .false.
    else
        i = 1
        do while (i <= n)
            call random_number(v1)
            call random_number(v2)
            v1 = 2.0_pr * v1 - 1.0_pr
            v2 = 2.0_pr * v2 - 1.0_pr
            rsq = v1*v1 + v2*v2

            if (rsq > 0.0_pr .and. rsq < 1.0_pr) then
                rsq = sqrt(-2.0_pr * log(rsq) / rsq)
                harvest(i) = v1 * rsq
                if (i < n) then
                    g(i+1) = v2 * rsq
                end if
                i = i + 2
            end if
        end do
        gaus_stored = .true.
    end if

end subroutine gasdev_v

!##################################################################################################
!     Not used / Not implemented
!##################################################################################################


    !subroutine create_file_name(prefix, num, suffix, filename)
    !    character(len=8)                :: fmt  ! Format descriptor
    !    character(len=12)               :: x1   ! Temporary string for formatted real
    !    character(len=*), intent(in)    :: prefix, suffix
    !    character(len=*)                :: filename
    !    real(kind=pr), intent(in)       :: num  ! Input real number
    !
    !    fmt = '(I10)'  ! Format integer
    !    write(x1, fmt) num  ! Convert real to string
    !
    !    ! Trim spaces in formatted number
    !    x1 = adjustl(trim(x1))
    !
    !    ! Concatenate strings
    !    filename = prefix // trim(x1) // suffix
    !
    !end subroutine create_file_name

    !subroutine Coulomb(particle_distance_squared, particle_separation,  force_contribution, E_potential, pressure_virial &
    !    , potential_cutoff)
    !    real(pr), intent(in)       :: particle_distance_squared, particle_separation(3)
    !    real(pr), intent(out)      :: force_contribution(3)
    !    real(pr), intent(inout)    :: E_potential, pressure_virial
    !    real(pr), intent(in)       :: potential_cutoff
    !
    !
    !end subroutine Coulomb

END MODULE