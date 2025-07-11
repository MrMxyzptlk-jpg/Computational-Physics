MODULE MD_rescaleMod
    use variablesMod
    use propertiesMod
    implicit none

CONTAINS

subroutine thermostat_rescale()
    real(pr)    :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities) / (3.0_pr * real(num_atoms,pr))
    scaling_factor = sqrt(ref_Temp / instant_Temp)

    velocities = velocities*scaling_factor

end subroutine thermostat_rescale

END MODULE MD_rescaleMod