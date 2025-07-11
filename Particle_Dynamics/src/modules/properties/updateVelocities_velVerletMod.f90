MODULE updateVelocities_velVerletMod
    use precisionMod
    use propertiesMod
    use variablesMod

CONTAINS

subroutine update_velocities_velVer()

    velocities = velocities + dt * 0.5_pr*(previous_forces + forces)

end subroutine update_velocities_velVer

END MODULE updateVelocities_velVerletMod