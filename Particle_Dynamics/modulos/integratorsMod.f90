MODULE integratorsMod
    use precisionMod
    use subroutinesMod
    use constantsMod
    use mzranMod
    use parsingMod
    use writing2filesMod
    use updatePositionsMod
    use observablesMod
    use initializationsMod
    implicit none
contains

subroutine get_measurements(i_measure)      ! Should be in observables module
    integer(int_large)      :: i_measure

    call get_observables(velocities, Energies(2,i_measure), pressures(i_measure), temperatures(i_measure))
    if ((.not.transitory) .and. do_mean_sqr_displacement) call update_msd(positions, meanSqrDisplacement)
    if (do_structure_factor) then
        call get_structure_factor(positions, structure_factor(i_measure), reciprocal_vec)
        call write_tasks(real(i_measure*measuring_jump,pr)*dt, positions, velocities, energies(:,i_measure) &
            , pressures(i_measure), temperatures(i_measure), structure_factor(i_measure))
    else
        call write_tasks(real(i_measure*measuring_jump,pr)*dt, positions, velocities, energies(:,i_measure) &
            , pressures(i_measure), temperatures(i_measure), structure_factor(1))
    end if

end subroutine get_measurements

subroutine velVerlet_step(i_measure)
    integer(int_large)                  :: i_measure

    call update_positions_velVer(positions, velocities, forces)
    previous_forces = forces
    if (do_linkCell) call create_links(positions)
    call get_forces(positions, forces, Energies(1,i_measure), pressures(i_measure), pair_corr)  ! If measure = .False. the observables are ignored

    call update_velocities_velVer(velocities, forces, previous_forces)

end subroutine velVerlet_step

subroutine Brownian_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    if (do_linkCell) call create_links(positions)
    call get_forces(positions, forces, Energies(1, i_measure), pressures(i_measure), pair_corr)

    call update_positions_Brownian(positions, forces)

end subroutine Brownian_step

subroutine MC_step(i_measure)
    integer(int_large)      :: i_measure

    if (do_linkCell) call create_links(positions)
    call get_forces(positions, forces, Energies(1,i_measure), pressures(i_measure), pair_corr)

    call update_positions_MC(positions, Energies(1,i_measure), MC_accepted)

end subroutine MC_step

END MODULE integratorsMod