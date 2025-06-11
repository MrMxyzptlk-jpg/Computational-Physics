MODULE initializations
    use precision
    use subrutinas
    use parsing
    use linkedLists
    implicit none

    real(pr), dimension(:), allocatable     :: Pressures, Temperatures, pair_corr, meanSqrDisplacement
    real(pr), dimension(:,:), allocatable   :: positions, velocities, forces, Energies
    real(pr)                                :: reciprocal_vec(3), structure_factor
    integer(int_large)                      :: transitory_minIndex, MC_accepted
    logical                                 :: do_linkCell = .False.

    abstract interface
        subroutine init_pos()
        end subroutine init_pos

        subroutine force_sub(positions, forces, E_potential, pressure_virial, pair_corr)
            use precision
            real(pr), intent(in)    :: positions(:,:)
            real(pr), intent(out)   :: forces(:,:)
            real(pr), intent(out)   :: E_potential, pressure_virial
            real(pr), intent(inout) :: pair_corr(:)
        end subroutine force_sub

        subroutine pairCorr_sub(positions, pair_corr)
            use precision
            real(pr), intent(in)        :: positions(:,:)
            real(pr), intent(inout)     :: pair_corr(:)
        end subroutine pairCorr_sub

        subroutine thermo(velocities)
            use precision
            real(pr), allocatable, intent(inout)    :: velocities(:,:)
        end subroutine thermo
    end interface

    procedure(init_pos), pointer         :: initialize_positions    => null()
    procedure(force_sub), pointer    :: get_forces              => null()
    procedure(pairCorr_sub), pointer :: get_pair_correlation    => null()
    procedure(thermo), pointer       :: thermostat_chosen       => null()

CONTAINS

subroutine init_structure()

    select case (structure)
        case ("FCC")
            initialize_positions =>  initialize_positions_FCC
            if (density > 0._pr) lattice_constant = (4._pr/density)**(1._pr/3._pr)
        case ("BCC")
            initialize_positions =>  initialize_positions_BCC
            if (density > 0._pr) lattice_constant = (2._pr/density)**(1._pr/3._pr)
        case ("random")
            initialize_positions =>  initialize_positions_random
            if (density > 0._pr) lattice_constant = (1._pr/density)**(1._pr/3._pr)
            if (density > 0._pr) then
                print*, "Random structure selected -> Chosen density is being ignored. Using density = num_atoms/vol instead"
            end if
        case default
            initialize_positions =>  initialize_positions_FCC
            if (density > 0._pr) lattice_constant = (4._pr/density)**(1._pr/3._pr)
    end select

end subroutine init_structure

subroutine init_potential()

    select case (type)
        !case ("Coulomb")
        !    potential =>  Coulomb  ! Not implemented yet
        case ("lannard_jones")
            potential =>  Lennard_Jones
            potential_function => Lennard_Jones_potential
        case default
            potential =>  Lennard_Jones
            potential_function => Lennard_Jones_potential
    end select

end subroutine init_potential

subroutine init_thermostat()

    select case (thermostat_type)
        case ("rescale")
            thermostat_chosen =>  thermostat_rescale
        case ("Berendsen")
            thermostat_chosen =>  thermostat_Berendsen
        case default
            thermostat_chosen =>  thermostat_rescale
    end select

end subroutine init_thermostat

subroutine init_observables()

    select case(integrator)
        case ('velocity-Verlet')
            if (save_transitory) then
                transitory_minIndex = -int(transitory_steps/measuring_jump)
                allocate(Energies(2,transitory_minIndex:measuring_steps)) ! Energies = (E_potential, E_kinetic)
                allocate(Pressures(transitory_minIndex:measuring_steps), Temperatures(transitory_minIndex:measuring_steps))
            else
                allocate(Energies(2,0:measuring_steps)) ! Energies = (E_potential, E_kinetic)
                allocate(Pressures(0:measuring_steps), Temperatures(0:measuring_steps))
            end if
            allocate(velocities(size(positions,1),size(positions,2)))
            allocate(forces(size(positions,1),size(positions,2)))
        case('Monte-Carlo')
            if (save_transitory) then
                transitory_minIndex = -int(transitory_steps/measuring_jump)
                allocate(Energies(1,transitory_minIndex:measuring_steps)) ! Energies = (E_potential, E_kinetic)
            else
                allocate(Energies(1,0:measuring_steps)) ! Energies = (E_potential, E_kinetic)
            end if
            allocate(Pressures(1), Temperatures(1))
    end select

end subroutine init_observables

subroutine init_summation()

    if (summation == "linked-lists") call check_linkCell(do_linkCell)

    if (do_linkCell) then
        get_forces =>  get_forces_linkedList
        call create_maps()
        call create_links(positions)
        if (integrator == 'Monte-Carlo' .and. do_pair_correlation) get_pair_correlation => get_pair_correlation_linkedlist
    else
        summation = "all-vs-all"
        get_forces =>  get_forces_allVSall
        if (integrator == 'Monte-Carlo' .and. do_pair_correlation) get_pair_correlation => get_pair_correlation_allVSall
    end if

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

subroutine initialize_positions_random() ! Not debugged
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

end subroutine initialize_positions_random

subroutine initialize_positions_FCC()
    integer(int_large)                 :: supercell_atoms
    integer(int_large)                 :: atom_id
    integer(int_medium)                :: h, k, l, b
    integer(int_small)                 :: i

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

end subroutine initialize_positions_FCC

subroutine initialize_positions_BCC()
    integer(int_large)                 :: supercell_atoms
    integer(int_large)                 :: atom_id
    integer(int_medium)                :: h, k, l, b
    integer(int_small)                 :: i

    ! BCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,2) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.5_pr  &  ! Atom 2
    ], [3,2])

    supercell_atoms = 2
    do i = 1_int_small, size(cell_dim)
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

end subroutine initialize_positions_BCC

subroutine initialize_velocities_random()
    real(pr)                                :: velocity_average(size(velocities,1))
    integer                                 :: i

    velocities = reshape( [ (rmzran() - 0.5_pr, i = 1, size(velocities)) ], shape(velocities) )

    velocity_average = sum(velocities,2)/real(num_atoms,pr)

    do i = 1, 3
        velocities(i,:) = (velocities(i,:) - velocity_average(i))   ! Removing center of mass displacement
    end do

end subroutine initialize_velocities_random

subroutine initialize_velocities_Maxwell()
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
    velocity_average = sum(velocities, dim=2) / real(num_atoms, pr)

    do i = 1, 3
        velocities(i,:) = velocities(i,:) - velocity_average(i)
    end do

    deallocate(tmp)

end subroutine initialize_velocities_Maxwell

subroutine initialize_velocities()

    select case (initial_velocities)
        case('random')
            call initialize_velocities_random()
        case('Maxwell')
            call initialize_velocities_Maxwell()
    end select

end subroutine initialize_velocities

END MODULE initializations