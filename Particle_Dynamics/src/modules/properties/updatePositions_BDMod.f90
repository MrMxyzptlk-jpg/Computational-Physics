MODULE updatePositions_BDMod
    use variablesMod
    use potentialsMod
    use subroutinesMod
    use propertiesMod
    use randomMod
    implicit none

CONTAINS

subroutine update_positions_Brownian()
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

END MODULE updatePositions_BDMod