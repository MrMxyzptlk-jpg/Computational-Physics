MODULE phase_runMod
    use precisionMod
    use variablesMod
    use propertiesMod
    use measurementsMod
    use writing2filesMod
    use integratorsMod, only : integrator_step
    implicit none

    integer(int_large)  :: i_measure

CONTAINS

subroutine phase_run()

    if (save_transitory) then
        i_measure = transitory_minIndex
        measure = .true.
    else
        i_measure = 0
    end if

    print*, "Integrator: ", integrator

    call open_files(reciprocal_vec)
        call get_forces(Energies(1,i_measure), pressures(i_measure), pair_corr)

        if (save_transitory) call get_measurements(i_measure)

        call run_transitory()

        transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
        measure     = .False.
        i_measure   = 0

        call run_definitive()

    call close_files()

end subroutine phase_run

subroutine run_transitory()
    integer(int_large)  :: i, j

    do i = -transitory_steps/thermostat_steps , -1, 1
        do j = 1, thermostat_steps
            if (save_transitory) call check_measuring(i*thermostat_steps + j, i_measure)    ! Checks if there will be measurements in this iteration

            call integrator_step(i_measure)

            if (measure) call get_measurements(i_measure)

        end do
        call thermostat_chosen()
    end do

end subroutine run_transitory

subroutine run_definitive()
    integer(int_large)  :: i

    do i = 1 , real_steps
        call check_measuring(i, i_measure)    ! Checks if there will be measurements in this iteration

        call integrator_step(i_measure)

        if (measure) call get_measurements(i_measure)

        if ((ensemble=='NVT').and.(mod(i,thermostat_steps)==0)) call thermostat_chosen()
    end do

end subroutine run_definitive

END MODULE phase_runMod