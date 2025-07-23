MODULE velVerlet_stepMod
    use precisionMod
    use subroutinesMod
    use parsingMod
    use updatePositionsMod
    use propertiesMod
    use forcesMod, only : get_forces, create_links
    use updateVelocities_velVerletMod
    implicit none

contains

subroutine velVerlet_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    call update_positions_velVer()
    previous_forces = forces
    if (do_linkCell) call create_links()
    call get_forces(Energies(1,i_measure), pressures(i_measure), pair_corr)  ! If measure = .False. the observables are ignored

    call update_velocities_velVer()

end subroutine velVerlet_step

END MODULE velVerlet_stepMod