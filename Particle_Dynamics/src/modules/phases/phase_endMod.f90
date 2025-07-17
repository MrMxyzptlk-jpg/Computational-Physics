MODULE phase_endMod
    use precisionMod
    use parsingMod
    use checkParsingMod
    use omp_lib
    use writing2filesMod
    use initializationsMod
    use integratorsMod
    use propertiesMod
    use observablesMod
    use measurementsMod
    implicit none

CONTAINS

subroutine phase_end(CPU_elapsed_time)
    real(pr), intent(in)    :: CPU_elapsed_time

    call open_observables_file(reciprocal_vec)

        if (save_transitory) call end_transitory()

        if (save_observables) call end_definitive()
    call close_observables_file()

    if (do_pair_correlation) then
        call normalize_pair_correlation(pair_corr)
        call write_pair_corr(pair_corr)
    end if

    if (do_mean_sqr_displacement) then
        call normalize_msd(meanSqrDisplacement)
        call write_msd(meanSqrDisplacement)
    end if

    if (save_state) call write_stateXML()

    call write_output(CPU_elapsed_time, energies(:,1:)/real(num_atoms,pr), pressures(1:), temperatures(1:), structure_factor(1:))

end subroutine phase_end


subroutine end_transitory()
    integer(int_large)  :: i, j, i_measure

    measure     = .True.
    i_measure = transitory_minIndex

    do i = -transitory_steps/thermostat_steps , -1, 1
        do j = 1, thermostat_steps
            if (save_transitory) call check_measuring(i*thermostat_steps + j, i_measure)    ! Checks if there will be measurements in this iteration

            if (do_structure_factor) then
                call write_tasks(real(i_measure*measuring_jump,pr)*dt, energies(:,i_measure)/real(num_atoms,pr) &
                    , pressures(i_measure), temperatures(i_measure), structure_factor(i_measure))
            else
                call write_tasks(real(i_measure*measuring_jump,pr)*dt, energies(:,i_measure)/real(num_atoms,pr) &
                    , pressures(i_measure), temperatures(i_measure), structure_factor(1))
            end if

        end do
    end do

end subroutine end_transitory

subroutine end_definitive()
    integer(int_large)  :: i, i_measure

    measure     = .True.
    i_measure   = 0

    do i = 1 , real_steps
        call check_measuring(i, i_measure)    ! Checks if there will be measurements in this iteration

        if (do_structure_factor) then
            call write_tasks(real(i_measure*measuring_jump,pr)*dt, energies(:,i_measure)/real(num_atoms,pr) &
                , pressures(i_measure), temperatures(i_measure), structure_factor(i_measure))
        else
            call write_tasks(real(i_measure*measuring_jump,pr)*dt, energies(:,i_measure)/real(num_atoms,pr) &
                , pressures(i_measure), temperatures(i_measure), structure_factor(1))
        end if
    end do

end subroutine end_definitive

END MODULE phase_endMod