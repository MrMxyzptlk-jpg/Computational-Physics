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

    call open_XYZ_file()
        call get_forces(Energies(1,i_measure), pressures(i_measure), pair_corr)

        if (save_transitory) then
            call get_measurements(i_measure)
            if (save_positions) call write_XYZfile(real(i_measure*measuring_jump,pr)*dt)
        end if

        print*, "Starting transitory run"
        call run_transitory()

        transitory  = .False. ! Flag to avoid calculations and saving variables during the transitory steps. False means the calculations are now NOT transitory
        measure     = .False.
        i_measure   = 0

        print*, "Starting definitive run"
        call run_definitive()

        if (debugg) call print_debuggs()

    call close_XYZ_file()

end subroutine phase_run

subroutine run_transitory()
    integer(int_large)  :: i, j

    do i = -transitory_steps/thermostat_steps , -1, 1
        do j = 1, thermostat_steps
            if (save_transitory) call check_measuring(i*thermostat_steps + j, i_measure)    ! Checks if there will be measurements in this iteration

!            if (interactions == "Coulomb" .and. integrator == "Monte-Carlo") call update_reciprocalCharges()   ! This is not an often enough recalculation of the reciprocal charges. Too much rounding error

            call integrator_step(i_measure)

            if (measure) then
                call get_measurements(i_measure)
                if (save_positions) call write_XYZfile(real(i_measure*measuring_jump,pr)*dt)
            end if

        end do
        call thermostat_chosen()
    end do

end subroutine run_transitory

subroutine run_definitive()
    integer(int_large)  :: i

    do i = 1 , real_steps
        call check_measuring(i, i_measure)    ! Checks if there will be measurements in this iteration

!        if (interactions == "Coulomb" .and. integrator == "Monte-Carlo") call update_reciprocalCharges()   ! This is not an often enough recalculation of the reciprocal charges. Too much rounding error

        call integrator_step(i_measure)

        if (measure) then
            call get_measurements(i_measure)
            if (save_positions) call write_XYZfile(real(i_measure*measuring_jump,pr)*dt)
        end if

        if ((ensemble=='NVT').and.(mod(i,thermostat_steps)==0)) call thermostat_chosen()
    end do

end subroutine run_definitive

subroutine print_debuggs()

        if (interactions == "Coulomb") then
            print*, "σ =", sigma
            print*, "Potential Contribution (Real)          :", E_potential_real/real(measuring_steps*num_atoms,pr)
            print*, "Potential Contribution (Reciprocal)    :", E_potential_reciprocal/real(measuring_steps*num_atoms,pr)
            print*, "Potential Contribution (Self Term)     :", Ewald_selfTerm/real(num_atoms,pr)
            print*, "Potential Contribution (Jelium Term)   :", Ewald_jeliumTerm/real(num_atoms,pr)
            print*, "Potential Contribution (Total)         :"&
                , ((E_potential_real + E_potential_reciprocal)/real(measuring_steps,pr) - Ewald_selfTerm - Ewald_jeliumTerm) &
                /real(num_atoms,pr)
            print*, "Net charge =", sum(charges(:))
        end if

end subroutine print_debuggs

END MODULE phase_runMod