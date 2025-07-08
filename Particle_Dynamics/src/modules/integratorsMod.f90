MODULE integratorsMod
    use precisionMod
    use subroutinesMod
    use parsingMod
    use updatePositionsMod
    use propertiesMod
    use forcesMod, only : get_forces, create_links
    implicit none

    abstract interface
        subroutine integrate (step)
            use precisionMod
            integer(int_large), intent(in) :: step
        end subroutine integrate
    end interface
    procedure(integrate), pointer   :: integrator_step  => null()

contains

subroutine velVerlet_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    call update_positions_velVer()
    previous_forces = forces
    if (do_linkCell) call create_links(positions)
    call get_forces(Energies(1,i_measure), pressures(i_measure), pair_corr)  ! If measure = .False. the observables are ignored

    call update_velocities_velVer()

end subroutine velVerlet_step

subroutine Brownian_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    if (do_linkCell) call create_links(positions)
    call get_forces(Energies(1, i_measure), pressures(i_measure), pair_corr)

    call update_positions_Brownian()

end subroutine Brownian_step

subroutine MC_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    if (do_linkCell) call create_links(positions)
    call get_forces(Energies(1,i_measure), pressures(i_measure), pair_corr)

    call update_positions_MC(Energies(1,i_measure), MC_accepted)

end subroutine MC_step

END MODULE integratorsMod