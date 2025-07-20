MODULE initializationsMod
    use precisionMod
    use constantsMod
    use randomMod
    use variablesMod
    use potentialsMod
    use forcesMod
    use thermostatsMod
    use subroutinesMod
    use parsingMod
    use observablesMod
    use propertiesMod
    use integratorsMod
    implicit none

    abstract interface
        subroutine init_pos()
        end subroutine init_pos

        subroutine pairCorr_sub(positions, pair_corr)
            use precisionMod
            real(pr), intent(in)        :: positions(:,:)
            real(pr), intent(inout)     :: pair_corr(:)
        end subroutine pairCorr_sub
    end interface

    procedure(init_pos), pointer     :: init_positions      => null()

CONTAINS

subroutine init_structure()

    select case (structure)
        case ("FCC")
            init_positions =>  init_positions_FCC
            lattice_constant = (4._pr/density)**(1._pr/3._pr)
        case ("BCC")
            init_positions =>  init_positions_BCC
            lattice_constant = (2._pr/density)**(1._pr/3._pr)
        case ("random")
            init_positions =>  init_positions_random
            cell_dim = 1_int_medium
            lattice_constant = (1._pr/density)**(1._pr/3._pr)
    end select

    if (state == 'fromFile') init_positions =>  init_positions_fromFile

end subroutine init_structure

subroutine init_variables()

    ! Define conversion factors to adimensionalize the variables
    conversion_factors(1) = sigma                   ! Distance      [Bohr = a₀]
    conversion_factors(2) = sigma*sqrt(mass/delta)  ! Time          [ℏ/Eh = tₐ]
    conversion_factors(3) = 1._pr                   ! Temperature   [Eh / kB]
    conversion_factors(4) = delta                   ! Energy        [hartree = Eh]
    conversion_factors(5) = sqrt(delta/mass)        ! Velocity      [a₀ / tₐ]
    conversion_factors(6) = delta/(sigma**3)        ! Pressure      [hartree / bohr³]

    density = density*conversion_factors(1)**3
    lattice_constant = lattice_constant/conversion_factors(1)
    periodicity = cell_dim*lattice_constant
    volume = product(periodicity)

    print*, "Periodicity = ", periodicity*conversion_factors(1)
    print*, "Volume", volume*conversion_factors(1)**3

    ref_Temp = ref_Temp*conversion_factors(3)    ! Already non-dimensional!!
    radius_cutoff = radius_cutoff/conversion_factors(1)
    dt = dt/conversion_factors(2)

    dtdt = dt*dt
    radius_cutoff_squared = radius_cutoff*radius_cutoff

    if (interactions == "lennard-jones") then
        potential_cutoff = 0._pr    ! Must be initialized as 0 because it will be used in the potential_function()
        potential_cutoff = potential_function(0, 0, radius_cutoff_squared) ! the first two arguments are there in case we want to implement non-equivalent particles latter on
    end if

    transitory = .True.    ! Flag to avoid calculations and saving variables during the transitory steps

    if(do_pair_correlation) then
        pair_corr_cutoff = pair_corr_cutoff/conversion_factors(1)
        pair_corr_cutoff_sqr = pair_corr_cutoff*pair_corr_cutoff
        dr = pair_corr_cutoff/real(pair_corr_bins,pr)
    end if

    call mzran_init() ! Initialize random number generator

    measure = .False.
    measuring_steps = real_steps/measuring_jump

end subroutine init_variables

subroutine init_potential()

    select case (interactions)
        case ("lennard-jones")
            potential =>  Lennard_Jones
            potential_function => Lennard_Jones_potential
        case ("Coulomb")
            sigma_sqr  = sigma*sigma
            potential =>  Coulomb_Ewald_realSpace
            potential_function =>  Coulomb_realSpace
            potential_function_reciprocal =>  Coulomb_reciprocalSpace
        case ("reaction_field")
            potential =>  reaction_field
            potential_function => reaction_field_potential
    end select

end subroutine init_potential

subroutine init_thermostat()

    if (integrator=="velocity-Verlet") then
        select case (thermostat_type)
            case ("rescale")
                thermostat_chosen =>  thermostat_rescale
            case ("Berendsen")
                thermostat_chosen =>  thermostat_Berendsen
            case default
                thermostat_chosen =>  thermostat_rescale
        end select
    end if

    if (integrator=="Monte-Carlo") thermostat_chosen =>  update_random_step

end subroutine init_thermostat

subroutine init_observables()
    if (save_transitory) then
        transitory_minIndex = -int(transitory_steps/measuring_jump)
        allocate(Energies(2,transitory_minIndex:measuring_steps)) ! Energies = (E_potential, E_kinetic)
        allocate(Pressures(transitory_minIndex:measuring_steps), Temperatures(transitory_minIndex:measuring_steps))
        if (do_structure_factor) allocate(structure_factor(transitory_minIndex:measuring_steps))
    else
        allocate(Energies(2,0:measuring_steps)) ! Energies = (E_potential, E_kinetic)
        allocate(Pressures(0:measuring_steps), Temperatures(0:measuring_steps))
        if (do_structure_factor) allocate(structure_factor(0:measuring_steps))
    end if
    allocate(forces(3,num_atoms))

    if (integrator == 'velocity-Verlet') then
        allocate(previous_forces(3,num_atoms))
        if (state == 'fromScratch') allocate(velocities(3,num_atoms))
    end if
    if (.not. do_structure_factor) allocate(structure_factor(1))

end subroutine init_observables

subroutine init_integrator()

    select case(integrator)
        case("velocity-Verlet")
            integrator_step => velVerlet_step
        case("Monte-Carlo")
            integrator_step => MC_step
            if (summation == "Ewald") then
                get_E_potential_contribution  => get_E_potential_contribution_Ewald
                update_potential_contribution => update_potential_contribution_Ewald
            else
                get_E_potential_contribution => get_E_potential_contribution_normal
                update_potential_contribution => update_potential_contribution_normal
            end if
        case("Brownian")
            integrator_step => Brownian_step
    end select

end subroutine init_integrator

subroutine init_summation()

    select case (summation)
        case ("linked-lists")
            call check_linkCell(do_linkCell)
            if (do_linkCell) then
                get_forces =>  get_forces_linkedList
                call create_maps()
                call create_links()
            else
                summation = "all-vs-all"
                get_forces =>  get_forces_allVSall
            end if
        case ("all-vs-all")
            get_forces =>  get_forces_allVSall
        case ("Ewald")
            get_forces =>  get_forces_Ewald
    end select

end subroutine init_summation

subroutine init_tasks()

    if (do_structure_factor) call get_reciprocal_vec(Miller_index, reciprocal_vec)

    if (do_mean_sqr_displacement) call initialize_msd(meanSqrDisplacement)

    if (do_pair_correlation) then
        allocate(pair_corr(pair_corr_bins))
        pair_corr = 0._pr
    else
        allocate(pair_corr(1)) ! Must be allocated to avoid issues in parallelized subroutines
        pair_corr = 0._pr
    end if

end subroutine init_tasks

subroutine init_positions_random() ! Not debugged
    integer                             :: i, j, accepted, attempts
    real(pr)                            :: r2_min, r2, dx, dy, dz
    real(pr)                            :: r_min
    real(pr), allocatable               :: pos_trial(:,:)

    allocate(positions(3, num_atoms))
    positions = 0._pr
    r_min = sigma*0.5_pr   ! Minimum allowed distance between particles
    r2_min = r_min**2

    accepted = 0
    attempts = 0

    do while (accepted < num_atoms)
        ! Generate random position in the supercell
        pos_trial = reshape( [ (rmzran() * periodicity(i), i = 1, 3) ], [3, 1] )

        ! Check against all previously accepted positions
        do j = 1, accepted
            dx = pos_trial(1,1) - positions(1,j)
            dy = pos_trial(2,1) - positions(2,j)
            dz = pos_trial(3,1) - positions(3,j)

            r2 = dx*dx + dy*dy + dz*dz
            if (r2 < r2_min) exit  ! Too close → reject
        end do

        ! If we made it through the loop, the point is valid
        if (j > accepted) then
            accepted = accepted + 1
            positions(:,accepted) = pos_trial(:,1)
        end if

        attempts = attempts + 1
        if (attempts > 100000) then
            print *, "ERROR: Could not place all atoms with minimum spacing."
            stop
        end if
        print*, "Attempt:", attempts, " Accepted:", accepted
    end do

    print*, "Random positions initialized after ", attempts, "attempts"

end subroutine init_positions_random

subroutine init_positions_FCC()
    integer(int_large)      :: supercell_atoms
    integer(int_large)      :: atom_id
    integer                 :: i, h, k, l, b

    ! FCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,4) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.0_pr, &  ! Atom 2
        0.5_pr, 0.0_pr, 0.5_pr, &  ! Atom 3
        0.0_pr, 0.5_pr, 0.5_pr   & ! Atom 4
    ], [3,4])

    supercell_atoms = 4
    do i = 1_int_small, size(cell_dim)
        supercell_atoms = supercell_atoms * cell_dim(i)
    end do
    if (supercell_atoms/=num_atoms .and. num_atoms/=0) then
        print'(a,I5,a,I5)', "Number of atoms mismatch: selected number of atoms = ", num_atoms &
            , " and atoms in the FCC supercell = ",supercell_atoms
        print*, "Ignoring number of atoms specified"
        num_atoms = supercell_atoms
    end if
    allocate(positions(3, num_atoms))

    atom_id = 1
    do h = 0_int_medium, cell_dim(1)-1_int_medium
        do k = 0_int_medium, cell_dim(2)-1_int_medium
            do l = 0_int_medium, cell_dim(3)-1_int_medium
                do b = 1, 4
                    positions(:, atom_id) = [real(h, pr), real(k, pr), real(l, pr)] + basis(:, b)
                    atom_id = atom_id + 1
                end do
            end do
        end do
    end do

    positions = positions*lattice_constant

end subroutine init_positions_FCC

subroutine init_positions_BCC()
    integer(int_large)      :: supercell_atoms
    integer(int_large)      :: atom_id
    integer                 :: i, h, k, l, b

    ! BCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,2) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.5_pr  &  ! Atom 2
    ], [3,2])

    supercell_atoms = 2
    do i = 1, size(cell_dim)
        supercell_atoms = supercell_atoms * cell_dim(i)
    end do
    if (supercell_atoms/=num_atoms .and. num_atoms/=0) then
        print'(a,I5,a,I5)', "Number of atoms mismatch: selected number of atoms = ", num_atoms &
            , " and atoms in the BCC supercell = ",supercell_atoms
        print*, "Ignoring number of atoms specified"
        num_atoms = supercell_atoms
    end if
    allocate(positions(3, num_atoms))

    atom_id = 1
    do h = 0_int_medium, cell_dim(1)-1_int_medium
        do k = 0_int_medium, cell_dim(2)-1_int_medium
            do l = 0_int_medium, cell_dim(3)-1_int_medium
                do b = 1, 2
                    positions(:, atom_id) = [real(h, pr), real(k, pr), real(l, pr)] + basis(:, b)
                    atom_id = atom_id + 1
                end do
            end do
        end do
    end do

end subroutine init_positions_BCC

subroutine init_positions_fromFile()
    integer     :: i, j
    real(pr)    :: parsed_periodicity(3)
    real(pr)    :: scaling_factor_vec(3), scaling_factor, scaling_tol

    call parse_stateXML(parsed_periodicity)   ! Getting the positions, num_atoms and periodicity from STATE.xml file

    ! Check consistent super-cell dimensions for re-scaling
    scaling_tol = 10._pr * epsilon(minval(periodicity) * conversion_factors(1))    ! Getting a tolerance base on the values of periodicity we are dealing with
    scaling_factor_vec = periodicity / parsed_periodicity   ! Note no re-dimensionalization of the periodicity to avoid another multiplication for the position re-scaling
    do i = 1, 3
        j = mod(i,3) + 1
        if (abs(scaling_factor_vec(i) - scaling_factor_vec(j)) > scaling_tol) then
            STOP "ERROR: Inconsistent periodicity rescaling from STATE.xml to input.xml"
        end if
    end do

    ! Re-scale positions to meet the new density from input.xml file. Note this does not propagate error after consecutive runs fromFile, because parameters remain unchanged
    scaling_factor = scaling_factor_vec(1)
    do i = 1, num_atoms
        positions(:,i) = positions(:,i) * scaling_factor
    end do

end subroutine init_positions_fromFile

subroutine init_velocities_random()
    real(pr)                                :: velocity_average(size(velocities,1))
    integer                                 :: i

    velocities = reshape( [ (rmzran() - 0.5_pr, i = 1, size(velocities)) ], shape(velocities) )

    velocity_average = sum(velocities,2)/real(num_atoms,pr)

    do i = 1, 3
        velocities(i,:) = (velocities(i,:) - velocity_average(i))   ! Removing center of mass displacement
    end do

end subroutine init_velocities_random

subroutine init_velocities_Maxwell()
    real(pr)                           :: velocity_average(3)
    real(pr), allocatable              :: tmp(:)
    integer                            :: i, n

    n = size(velocities)
    allocate(tmp(n))

    ! Generate Gaussian random values
    call gasdev_v(tmp)

    ! Convert to pr and reshape
    velocities = reshape([ (real(tmp(i), pr), i = 1, n) ], shape(velocities))

    ! Subtract center-of-mass velocity
    velocity_average = sum(velocities, 2) / real(num_atoms, pr)

    do i = 1, 3
        velocities(i,:) = velocities(i,:) - velocity_average(i)
    end do

    deallocate(tmp)

end subroutine init_velocities_Maxwell

subroutine init_velocities()
    select case (initial_velocities)
        case('Maxwell')
            call init_velocities_Maxwell()
        case('random')
            call init_velocities_random()
    end select
end subroutine init_velocities

subroutine init_internal_constants()

    if (interactions == "Coulomb") then
        if (.not. allocated(charges)) then; allocate(charges(num_atoms)); charges = 1._pr; end if
    end if

    if (integrator == 'Brownian') then
        if ((viscosity == 0._pr) .and. (reduced_viscosity == 0._pr)) then
            print*, "No viscosity of reduced viscosity specified. Setting reduced_viscosity = 1"
            reduced_viscosity = 1._pr
        else if (viscosity/=0._pr) then
            reduced_viscosity = 3._pr*pi*viscosity * sigma      ! 3πησ product
        else if(reduced_viscosity/=0._pr) then
            viscosity = 3._pr*pi * sigma / reduced_viscosity    ! 3πησ product
        end if

        diffusion_coeff = ref_Temp / reduced_viscosity          ! D₀ = k_B T / (3πησ)
        reduced_viscosity_inv = 1._pr / reduced_viscosity       ! 1 / (3πησ)
        brownian_stddev = sqrt(2.0_pr * diffusion_coeff * dt)
    end if

    Temp_factor = 2.0_pr / (3.0_pr * real(num_atoms,pr))
    Pressure_factor = 1._pr / (3._pr * volume)

end subroutine init_internal_constants

subroutine init_Ewald()
    integer     :: kx, ky, kz, kvec_count
    real(pr)    :: k_sqr
    logical     :: good_kvec

    ! Set all constants for Ewald summation
    k_periodicity = twoPi/periodicity
    Ewald_realFactor    = 2._pr  / (sigma * sqrt(pi))
    eightPi_over_volume = 8._pr*pi/volume
    twoPi_over_volume   = twoPi/volume
    radius_cutoff = 0.5_pr*minval(periodicity)  ! Set to half the minimum box width
    Ewald_selfTerm = sum(charges*charges)/ (sqrt(pi)*sigma)

    ! Get k-space factors for the force and potential energy contributions
    halfSigma_sqr = sigma*sigma / 4.0_pr

    ! Count valid vectors
    kvec_count = 0
    do kx = 0, kgrid(1)
        do ky = -kgrid(2), kgrid(2)
            do kz = -kgrid(3), kgrid(3)
                call check_kvec(kx, ky, kz, k_sqr, good_kvec)
                if (good_kvec) kvec_count = kvec_count + 1
            end do
        end do
    end do

    num_kvec = kvec_count
    allocate(k_vectors(num_kvec))
    kvec_count = 0

    ! Note that only 4 octants are considered and the others are accounted by symmetries later on.

    do kx = 0, kgrid(1)
        do ky = -kgrid(2), kgrid(2)
            do kz = -kgrid(3), kgrid(3)
                call check_kvec(kx, ky, kz, k_sqr, good_kvec)
                if (good_kvec) then
                    kvec_count = kvec_count + 1
                    k_vectors(kvec_count)%kx = kx
                    k_vectors(kvec_count)%ky = ky
                    k_vectors(kvec_count)%kz = kz
                    k_vectors(kvec_count)%k_sqr = k_sqr
                    k_vectors(kvec_count)%kvector = (/kx, ky, kz/) * k_periodicity

                    k_vectors(kvec_count)%kfactor = fourpi/volume * exp(-halfSigma_sqr * k_sqr) / k_sqr
                end if
            end do
        end do
    end do

end subroutine init_Ewald

subroutine init_reciprocalCharges()
    integer                 :: kx, ky, kz, i
    complex(pr)             :: eikx(0:kgrid(1), num_atoms)
    complex(pr)             :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr)             :: eikz(-kgrid(3):kgrid(3), num_atoms)

    allocate(reciprocal_charges(num_kvec))
    reciprocal_charges = 0._pr

    call get_all_expFactors(eikx, eiky, eikz)

    !$omp parallel private(kx, ky, kz, i) &
    !$omp shared(eikx, eiky, eikz, num_kvec, charges, k_vectors, reciprocal_charges) &
    !$omp default(none)

    !$omp do schedule(dynamic)
    do i = 1, num_kvec
        kx = k_vectors(i)%kx
        ky = k_vectors(i)%ky
        kz = k_vectors(i)%kz

        reciprocal_charges(i) = sum(charges(:)*eikx(kx,:)*eiky(ky,:)*eikz(kz,:))
    end do
    !$omp end do

    !$omp end parallel

end subroutine init_reciprocalCharges

!##################################################################################################
!     Not used / Not implemented
!##################################################################################################

!subroutine init_Ewald_old()
!    integer                         :: kx, ky, kz, kvec_count
!    real(pr)                        :: k_sqr, halfSigma_sqr, fourPi, k_periodicity(3)
!    type(kvector_data), allocatable :: temp_kvec(:)
!
!    halfSigma_sqr = sigma_sqr / 4.0_pr
!    fourPi = 4._pr * pi
!    k_periodicity = 2._pr*pi/periodicity
!
!
!    allocate(temp_kvec(product(2*kgrid+1)))  ! overestimate, shrink later
!    kvec_count = 0  ! Counter for the number of reciprocal lattice vectors in the first octant
!    do kx = 0, kgrid(1)
!        do ky = 0, kgrid(2)
!            do kz = 0, kgrid(3)
!                if (kx == 0 .and. ky == 0 .and. kz == 0) cycle
!                ! Only store kx > 0 or (kx==0 and ky>=0 and kz>=0) to avoid double-counting
!                if (kx > 0 .or. (kx==0 .and. ky >= 0 .and. kz >= 0)) then
!                    kvec_count = kvec_count + 1
!                    temp_kvec(kvec_count)%kvec = real((/kx, ky, kz/),pr)*k_periodicity
!                    k_sqr = sum((temp_kvec(kvec_count)%kvec)**2)
!                    temp_kvec(kvec_count)%k_squared = k_sqr
!                    temp_kvec(kvec_count)%k_factor = exp(-k_sqr*halfSigma_sqr) / k_sqr
!    !!                if (kx == 0 .and. ky == 0 .and. kz == 1) temp_kvec(kvec_count)%k_factor = temp_kvec(kvec_count)%k_factor*0.5_pr
!                end if
!            end do
!        end do
!    end do
!
!    num_kvec = kvec_count
!    allocate(kvectors(num_kvec))
!    kvectors(:) = temp_kvec(:num_kvec)
!    deallocate(temp_kvec)
!
!end subroutine init_Ewald_old

END MODULE initializationsMod