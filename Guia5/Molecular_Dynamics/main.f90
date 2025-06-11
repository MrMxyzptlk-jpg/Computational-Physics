
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
! Room for improvement:
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
    use initializations
    implicit none

    real(pr), dimension(:,:), allocatable   :: previous_forces
    real(pr)                                :: CPU_t_start, CPU_t_end, CPU_elapsed_time
!    character (len=:), allocatable          :: filename, prefix, file_root_word, suffix
    integer(int_large)                      :: i, j, i_measure


!###################################################################################################
!   Set default values and parse input file (parsing module)
!###################################################################################################

    call set_defaults()
    call parse_input()

!##################################################################################################
!      Necessary definitions, pointers, initializations and conversion factors (initializations module unless specified otherwise)
!##################################################################################################

    call init_structure()
    call init_potential()
    call init_thermostat()
    call init_observables()
    call init_summation()
    call init_tasks()
    call initialize_parameters()                ! subrutinas module
    call initialize_XYZ_data()                  ! writing2flies module
    call initialize_positions()
    if (integrator /= 'Monte-Carlo') then
        call initialize_velocities()
        call initialize_rest()                  ! subrutinas module
        call thermostat_rescale(velocities)     ! subrutinas module
    end if

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

            do i = -transitory_steps/thermostat_steps , -1, 1
                do j = 1, thermostat_steps
                    if (save_transitory) call check_measuring(i*thermostat_steps + j)    ! Checks if there will be measurements in this iteration
                    if (measure) i_measure = i_measure + 1

                    call update_positions_velVer(positions, velocities, forces)
                    previous_forces = forces
                    if (do_linkCell) call create_links(positions)
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

                if (measure) then
                    if (do_linkCell) call create_links(positions)
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
        call normalize_pair_correlation(pair_corr)
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