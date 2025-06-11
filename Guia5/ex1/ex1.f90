
!***************************************************************
! Program: ex1.f90
! Purpose:
!
! Description:
!
! Input: All parameters are in ATOMIC UNITS. The program adimensionalized the problem in order to calculate, and re-dimnesionalizes the values when writing to files.
!
!
! Output:
!
!
!
! Room for improvement: add check to the random initialization of positions in order to avoid collisions. Check cutoff_radius units (and relation to sigma).
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precision
    use funciones
    use subrutinas
    use constantes
    use mzranmod
    use parsing
    use linkedLists
    use omp_lib
    use writing2files
    implicit none

    real(pr), dimension(:), allocatable     :: Pressures, Temperatures, pair_corr, meanSqrDisplacement
    real(pr)                                :: reciprocal_vec(3), structure_factor
    real(pr), dimension(:,:), allocatable   :: positions, velocities, forces, previous_forces, Energies
    real(pr)                                :: CPU_t_start, CPU_t_end, CPU_elapsed_time
!    character (len=:), allocatable          :: filename, prefix, file_root_word, suffix
    integer(int_large)                      :: i, j, i_measure, transitory_minIndex, MC_accepted
    logical                                 :: do_linkCell = .False.

    abstract interface
        subroutine init_pos(positions)
            use precision
            real(pr), allocatable, intent(out) :: positions(:,:)
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

    procedure(init_pos), pointer    :: initialize_positions => null()
    procedure(force_sub), pointer   :: get_forces => null()
    procedure(pairCorr_sub), pointer   :: get_pair_correlation => null()
    procedure(thermo), pointer      :: thermostat_chosen => null()

!###################################################################################################
!   Set default values and parse input file
!###################################################################################################

    call set_defaults()
    call parse_input()

!##################################################################################################
!      Necessary definitions, pointers, initializations and conversion factors
!##################################################################################################

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

    select case (type)
        case ("Coulomb")
            potential =>  Coulomb
        case ("lannard_jones")
            potential =>  Lennard_Jones
            potential_function => Lennard_Jones_potential
        case default
            potential =>  Lennard_Jones
            potential_function => Lennard_Jones_potential
    end select

    select case (thermostat_type)
        case ("rescale")
            thermostat_chosen =>  thermostat_rescale
        case ("Berendsen")
            thermostat_chosen =>  thermostat_Berendsen
        case default
            thermostat_chosen =>  thermostat_rescale
    end select

    call initialize_parameters()
    call initialize_XYZ_data()
    call initialize_positions(positions)

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
            call initialize_velocities(velocities)
        case('Monte-Carlo')
            if (save_transitory) then
                transitory_minIndex = -int(transitory_steps/measuring_jump)
                allocate(Energies(1,transitory_minIndex:measuring_steps)) ! Energies = (E_potential, E_kinetic)
            else
                allocate(Energies(1,0:measuring_steps)) ! Energies = (E_potential, E_kinetic)
            end if
            allocate(Pressures(1), Temperatures(1))
    end select

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

    if (do_structure_factor) then
        call get_reciprocal_vec(Miller_index, reciprocal_vec)
    end if

    if (do_mean_sqr_displacement) then
        call initialize_msd(meanSqrDisplacement)
    end if

    if (do_pair_correlation) then
        allocate(pair_corr(pair_corr_bins))
        pair_corr = 0._pr
    else
        allocate(pair_corr(1)) ! Must be allocated to avoid issues in parallelized subroutines
        pair_corr = 0._pr
    end if

    call initialize_rest()
    if (integrator /= 'Monte-Carlo') call thermostat_rescale(velocities)

!##################################################################################################
!      Start of the calculations
!##################################################################################################

    if (save_transitory) then
        i_measure = transitory_minIndex
    else
        i_measure = 0
    end if

    CPU_t_start = omp_get_wtime()

    if (integrator == 'velocity-Verlet') then

        call open_files(reciprocal_vec)
            call get_forces(positions, forces,  Energies(1,i_measure), pressures(i_measure), pair_corr)

            if (save_transitory) then
                call get_observables(velocities, Energies(2,i_measure), pressures(i_measure) &
                    , temperatures(i_measure))
                if (do_structure_factor) call get_structure_factor(positions, structure_factor, reciprocal_vec)
                call write_tasks(i_measure*dt, positions, velocities, energies(:,i_measure) &
                    , pressures(i_measure), temperatures(i_measure), structure_factor)
            end if

            if (save_transitory .and. save_positions) call write_XYZfile(0._pr, positions, velocities)
            do i = -transitory_steps/thermostat_steps , -1, 1
                do j = 1, thermostat_steps
                    if (save_transitory) call check_measuring(i*thermostat_steps + j)    ! Checks if there will be measurements in this iteration
                    if (measure) i_measure = i_measure + 1

                    call update_positions_velVer(positions, velocities, forces)
                    if (do_linkCell) call create_links(positions)
                    previous_forces = forces
                    call get_forces(positions, forces, Energies(1,i_measure), pressures(i_measure), pair_corr)  ! If measure = .False. the observables are ignored

                    call update_velocities_velVer(velocities, forces, previous_forces)

                    if (measure) then
                        call get_observables(velocities, Energies(2,i_measure), pressures(i_measure), temperatures(i_measure))
                        if (do_structure_factor) call get_structure_factor(positions, structure_factor, reciprocal_vec)
                        call write_tasks(real(i_measure,pr)*dt, positions, velocities, energies(:,i_measure), pressures(i_measure) &
                            , temperatures(i_measure), structure_factor)
                    end if
                end do
                call thermostat_chosen(velocities)
            end do

            transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
            measure     = .False.
            i_measure   = 0
            do i = 1 , real_steps
                call check_measuring(i)    ! Checks if there will be measurements in this iteration
                if (measure) i_measure = i_measure + 1

                call update_positions_velVer(positions, velocities, forces)
                if (do_linkCell) call create_links(positions)
                previous_forces = forces
                call get_forces(positions, forces, Energies(1,i_measure), pressures(i_measure), pair_corr)
                call update_velocities_velVer(velocities, forces, previous_forces)
                call get_observables(velocities, Energies(2,i_measure), pressures(i_measure), temperatures(i_measure))

                if (measure) then
                    if (do_structure_factor) call get_structure_factor(positions, structure_factor, reciprocal_vec)
                    if (do_mean_sqr_displacement) call update_msd(positions, meanSqrDisplacement)
                    call write_tasks(real(i_measure,pr)*dt, positions, velocities, energies(:,i_measure), pressures(i_measure)  &
                        , temperatures(i_measure), structure_factor)
                end if
                if ((ensemble=='NVT').and.(mod(i,thermostat_steps)==0)) call thermostat_chosen(velocities)
            end do
        call close_files()
    end if


    if (integrator == 'Monte-Carlo') then

        call open_files(reciprocal_vec)
            if (save_transitory) then
                if (do_structure_factor) call get_structure_factor(positions, structure_factor, reciprocal_vec)
                call write_tasks(i_measure*dt, positions, velocities, energies(:,i_measure), pressures(1), temperatures(1) &
                    , structure_factor)
            end if
            do i = -transitory_steps/MC_adjust_step , -1, 1
                do j = 1, MC_adjust_step
                    if (save_transitory) call check_measuring(i*MC_adjust_step + j)    ! Checks if there will be measurements in this iteration
                    if (measure) i_measure = i_measure + 1

                    call update_positions_random(positions, Energies(1,i_measure), MC_accepted)

                    if (measure) then
                        if (do_linkCell) call create_links(positions)
                        if (do_structure_factor) call get_structure_factor(positions, structure_factor, reciprocal_vec)
                        call write_tasks(real(i_measure,pr)*dt, positions, velocities, energies(:,i_measure), pressures(1) &
                            , temperatures(1) , structure_factor)
                    end if
                end do
                call update_random_step(MC_accepted)
            end do

            transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
            measure     = .False.
            i_measure   = 0
            do i = 1 , real_steps
                call check_measuring(i)    ! Checks if there will be measurements in this iteration
                if (measure) i_measure = i_measure + 1

                call update_positions_random(positions, Energies(1,i_measure), MC_accepted)
                if (do_linkCell) call create_links(positions)

                if (measure) then
                    if (do_pair_correlation) call get_pair_correlation(positions, pair_corr)
                    if (do_structure_factor) call get_structure_factor(positions, structure_factor, reciprocal_vec)
                    if (do_mean_sqr_displacement) call update_msd(positions, meanSqrDisplacement)
                    call write_tasks(real(i_measure,pr)*dt, positions, velocities, energies(:,i_measure), pressures(1) &
                        , temperatures(1), structure_factor)
                end if
            end do
        call close_files()
    end if

    if (do_pair_correlation) then
        call normalize_pair_correlation(pair_corr, real_steps)
        call write_pair_corr(pair_corr)
    end if

    if (do_mean_sqr_displacement) then
        call normalize_msd(meanSqrDisplacement)
        call write_msd(meanSqrDisplacement)
    end if

    CPU_t_end = omp_get_wtime()
    CPU_elapsed_time = CPU_t_end - CPU_t_start

    call write_output(CPU_elapsed_time, energies, pressures, temperatures)

end program ex1