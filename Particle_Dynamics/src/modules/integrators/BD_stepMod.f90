MODULE BD_stepMod
    use precisionMod
    use subroutinesMod
    use parsingMod
    use updatePositionsMod
    use propertiesMod
    use forcesMod, only : get_forces, create_links
    implicit none

contains

subroutine Brownian_step(i_measure)
    integer(int_large), intent(in) :: i_measure

    if (do_linkCell) call create_links(positions)
    call get_forces(Energies(1, i_measure), pressures(i_measure), pair_corr)

    call update_positions_Brownian()

end subroutine Brownian_step

END MODULE BD_stepMod