MODULE updatePositions_MCMod
    use variablesMod
    use potentialsMod
    use propertiesMod
    use randomMod
    implicit none
CONTAINS

subroutine update_positions_MC(E_potential)
    real(pr), intent(inout)                 :: E_potential
    real(pr)                                :: proposed_position(3), E_potential_new, E_potential_old, dE
    integer                                 :: random_particle_id, i, j

    do j = 1, num_atoms
        ! Pick a random particle
        random_particle_id = int(rmzran()*num_atoms) + 1

        ! Propose a displacement
        proposed_position = positions(:,random_particle_id) + (/(MC_delta*(rmzran() - 0.5_pr), i = 1, 3)/)

        ! Apply periodic boundary conditions
        proposed_position = modulo(proposed_position, periodicity(:))

        ! Compute potential energy contribution (E_new - E_old)
        call get_E_potential_contribution(random_particle_id, proposed_position, dE)

        ! Metropolis criterion
        if (dE <= 0._pr ) then
            call update_potential_contribution(random_particle_id, proposed_position, E_potential, dE)
        else if ( rmzran() < exp(-dE / ref_Temp)) then
            call update_potential_contribution(random_particle_id, proposed_position, E_potential, dE)
        end if
    end do

end subroutine update_positions_MC

END MODULE updatePositions_MCMod