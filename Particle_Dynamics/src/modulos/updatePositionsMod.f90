MODULE updatePositionsMod
    use parametersMod
    use potentialsMod
    use subroutinesMod
    use randomMod
    implicit none
CONTAINS

subroutine update_positions_MC(positions, E_potential, N_accepted)
    real(pr), dimension(:,:), intent(inout) :: positions
    real(pr), intent(inout)                 :: E_potential
    integer, intent(inout)                  :: N_accepted
    real(pr)                                :: random_displacement(3)
    real(pr)                                :: old_position(3), E_potential_old, E_potential_new, dE
    integer                                 :: random_particle_id, i, j

    do j = 1, num_atoms
        ! Pick a random particle
        random_particle_id = int(rmzran()*num_atoms) + 1
        old_position = positions(:,random_particle_id)

        ! Compute potential energy contribution
        call get_E_potential_contribution(positions, random_particle_id, E_potential_old)

        ! Propose a displacement
        random_displacement = (/(MC_delta*(rmzran() - 0.5d0), i = 1, 3)/)
        positions(:,random_particle_id) = old_position + random_displacement

        ! Apply periodic boundary conditions
        positions(:,random_particle_id) = modulo(positions(:,random_particle_id), periodicity(:))

        ! Compute new potential energy contribution
        call get_E_potential_contribution(positions, random_particle_id, E_potential_new)

        dE = E_potential_new - E_potential_old

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

subroutine update_positions_velVer(positions, velocities, forces)
    real(pr), dimension(:,:), intent(inout)    :: positions, velocities, forces

    positions = positions + velocities*dt + forces*dtdt*0.5_pr

    ! Apply periodic boundary conditions
    !positions = mod(positions, spread(periodicity, dim=2, ncopies=size(positions,2)))
    positions(1,:) = modulo(positions(1,:), periodicity(1))
    positions(2,:) = modulo(positions(2,:), periodicity(2))
    positions(3,:) = modulo(positions(3,:), periodicity(3))

end subroutine update_positions_velVer

subroutine update_positions_Brownian(positions, forces)
    real(pr), dimension(:,:), intent(inout)     :: positions, forces
    real(pr)                                    :: random_noise(3)
    integer                                     :: i, j

    ! Update positions with force drift + random Gaussian noise
    do i = 1, num_atoms
        random_noise = (/(rnd_normal(brownian_stddev), j=1, 3)/)
        positions(:,i) = positions(:,i) + reduced_viscosity_inv * forces(:,i) * dt + random_noise
        !print'(*(E9.3,2x))',reduced_viscosity_inv * forces(:,i) * dt, random_noise, brownian_stddev
    end do

    ! Apply periodic boundary conditions
    !positions = mod(positions, spread(periodicity, dim=2, ncopies=size(positions,2)))
    positions(1,:) = modulo(positions(1,:), periodicity(1))
    positions(2,:) = modulo(positions(2,:), periodicity(2))
    positions(3,:) = modulo(positions(3,:), periodicity(3))

end subroutine update_positions_Brownian

END MODULE updatePositionsMod