MODULE MC_stepMod
    use precisionMod
    use subroutinesMod
    use parsingMod
    use updatePositionsMod
    use propertiesMod
    use forcesMod, only : get_forces, create_links
    implicit none

contains

subroutine MC_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    if (do_linkCell) call create_links()

    call update_positions_MC(Energies(1,i_measure))
    if (measure .and. save_observables) call get_forces(Energies(1,i_measure), pressures(i_measure), pair_corr)

end subroutine MC_step

END MODULE MC_stepMod