
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
! Room for improvement:  further modularization for comprehension and ease of use should be done in the 'subrutinas' modue. That is, to subdivide it into smaller modules in the 'modulos' folder
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precisionMod
    use parsingMod
    use omp_lib
    use writing2filesMod
    use initializationsMod
    use integratorsMod
    use observablesMod
    implicit none

    real(pr)                                :: CPU_t_start, CPU_t_end, CPU_elapsed_time
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
    call initialize_parameters()                ! parametersMod module
    call init_positions()
    call init_thermostat()
    call init_observables()
    call init_tasks()
    call init_summation()
    if (integrator == 'velocity-Verlet') call initialize_velocities()
    call init_internal_constants()
    if (integrator == 'velocity-Verlet') call thermostat_rescale(velocities)     ! thermostatsMod module
    call initialize_XYZ_data()                  ! writing2fliesMod module

!##################################################################################################
!      Start of the calculations
!##################################################################################################

    if (save_transitory) then
        i_measure = transitory_minIndex
        measure = .true.
    else
        i_measure = 0
    end if

    CPU_t_start = omp_get_wtime()

    if (integrator == 'velocity-Verlet') then

        call open_files(reciprocal_vec)
            call get_forces(positions, forces,  Energies(1,i_measure), pressures(i_measure), pair_corr)

            if (save_transitory) call get_measurements(i_measure)

            do i = -transitory_steps/thermostat_steps , -1, 1
                do j = 1, thermostat_steps
                    if (save_transitory) call check_measuring(i*thermostat_steps + j, i_measure)    ! Checks if there will be measurements in this iteration

                    call velVerlet_step(i_measure)

                    if (measure) call get_measurements(i_measure)

                end do
                call thermostat_chosen(velocities)
            end do

            transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
            measure     = .False.
            i_measure   = 0
            do i = 1 , real_steps
                call check_measuring(i, i_measure)    ! Checks if there will be measurements in this iteration

                call velVerlet_step(i_measure)

                if (measure) call get_measurements(i_measure)

                if ((ensemble=='NVT').and.(mod(i,thermostat_steps)==0)) call thermostat_chosen(velocities)
            end do
        call close_files()
    end if

    if (integrator == 'Brownian') then

        call open_files(reciprocal_vec)
            call get_forces(positions, forces,  Energies(1,i_measure), pressures(i_measure), pair_corr)

            if (save_transitory) call get_measurements(i_measure)

            do i = -transitory_steps/thermostat_steps , -1, 1
                do j = 1, thermostat_steps
                    if (save_transitory) call check_measuring(i*thermostat_steps + j, i_measure)    ! Checks if there will be measurements in this iteration

                    call Brownian_step(i_measure)

                    if (measure) call get_measurements(i_measure)

                end do
            end do

            transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
            measure     = .False.
            i_measure   = 0
            do i = 1 , real_steps
                call check_measuring(i, i_measure)    ! Checks if there will be measurements in this iteration

                call Brownian_step(i_measure)

                if (measure) call get_measurements(i_measure)

            end do
        call close_files()
    end if

    if (integrator == 'Monte-Carlo') then

        call open_files(reciprocal_vec)
            call get_forces(positions, forces,  Energies(1,i_measure), pressures(i_measure), pair_corr)
            if (save_transitory) call get_measurements(i_measure)

            do i = -transitory_steps/MC_adjust_step , -1, 1
                do j = 1, MC_adjust_step
                    if (save_transitory) call check_measuring(i*MC_adjust_step + j, i_measure)    ! Checks if there will be measurements in this iteration

                    call MC_step(i_measure)

                    if (measure) call get_measurements(i_measure)

                end do
                call update_random_step(MC_accepted)
            end do

            transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
            measure     = .False.
            i_measure   = 0
            do i = 1 , real_steps
                call check_measuring(i,i_measure)    ! Checks if there will be measurements in this iteration

                call MC_step(i_measure)

                if (measure) call get_measurements(i_measure)
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

    call write_output(CPU_elapsed_time, energies(:,1:), pressures(1:), temperatures(1:), structure_factor(1:))

end program ex1