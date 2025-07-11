MODULE updatePositions_velVerletMod
    use variablesMod
    use potentialsMod
    use propertiesMod
    implicit none

CONTAINS

subroutine update_positions_velVer()

    positions = positions + velocities*dt + forces*dtdt*0.5_pr

    ! Apply periodic boundary conditions
    !positions = mod(positions, spread(periodicity, dim=2, ncopies=size(positions,2))) ! A different way of applying PBC
    positions(1,:) = modulo(positions(1,:), periodicity(1))
    positions(2,:) = modulo(positions(2,:), periodicity(2))
    positions(3,:) = modulo(positions(3,:), periodicity(3))

end subroutine update_positions_velVer

END MODULE updatePositions_velVerletMod