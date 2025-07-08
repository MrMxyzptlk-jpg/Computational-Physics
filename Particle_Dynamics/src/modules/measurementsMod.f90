MODULE measurementsMod
    use precisionMod
    use subroutinesMod
    use writing2filesMod
    use propertiesMod
    use observablesMod
    use initializationsMod
    implicit none
CONTAINS

subroutine get_measurements(i_measure)
    integer(int_large)      :: i_measure

    call get_observables(velocities, Energies(2,i_measure), pressures(i_measure), temperatures(i_measure))
    if ((.not.transitory) .and. do_mean_sqr_displacement) call update_msd(positions, meanSqrDisplacement)
    if (do_structure_factor) then
        call get_structure_factor(positions, structure_factor(i_measure), reciprocal_vec)
        call write_tasks(real(i_measure*measuring_jump,pr)*dt, energies(:,i_measure) &
            , pressures(i_measure), temperatures(i_measure), structure_factor(i_measure))
    else
        call write_tasks(real(i_measure*measuring_jump,pr)*dt, energies(:,i_measure) &
            , pressures(i_measure), temperatures(i_measure), structure_factor(1))
    end if

end subroutine get_measurements

END MODULE measurementsMod