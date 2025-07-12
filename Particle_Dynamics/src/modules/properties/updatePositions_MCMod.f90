MODULE updatePositions_MCMod
    use variablesMod
    use potentialsMod
    use propertiesMod
    use randomMod
    implicit none
CONTAINS

subroutine update_positions_MC(E_potential, N_accepted)
    real(pr), intent(inout)                 :: E_potential
    integer, intent(inout)                  :: N_accepted
    real(pr)                                :: random_displacement(3)
    real(pr)                                :: old_position(3), E_potential_new, dE
    integer                                 :: random_particle_id, i, j

    do j = 1, num_atoms
        ! Pick a random particle
        random_particle_id = int(rmzran()*num_atoms) + 1
        old_position = positions(:,random_particle_id)

        ! Propose a displacement
        random_displacement = (/(MC_delta*(rmzran() - 0.5d0), i = 1, 3)/)
        positions(:,random_particle_id) = old_position + random_displacement

        ! Apply periodic boundary conditions
        positions(:,random_particle_id) = modulo(positions(:,random_particle_id), periodicity(:))

        ! Compute new potential energy contribution
        call get_E_potential_contribution(random_particle_id, E_potential_new)

        dE = E_potential_new - previous_E_potential(j)

        ! Metropolis criterion
        if (dE <= 0._pr ) then
            if (measure .and. save_observables) E_potential = E_potential_new
            N_accepted = N_accepted + 1
            return
        else if ( rmzran() < exp(-dE / ref_Temp)) then
            if (measure .and. save_observables) E_potential = E_potential_new
            N_accepted = N_accepted + 1
            return
        else
            ! Revert move if trial not accepted
            positions(:,random_particle_id) = old_position
        end if
    end do

end subroutine update_positions_MC

END MODULE updatePositions_MCMod