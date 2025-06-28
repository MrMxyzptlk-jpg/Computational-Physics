MODULE thermostatsMod
    use parametersMod
    implicit none
CONTAINS

subroutine thermostat_rescale(velocities)
    real(pr), allocatable, intent(inout)    :: velocities(:,:)
    real(pr)                                :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities)/(3.0_pr * real(num_atoms,pr))
    scaling_factor = sqrt( ref_Temp / instant_Temp )

    velocities = velocities*scaling_factor

end subroutine thermostat_rescale

subroutine thermostat_Berendsen(velocities)
    real(pr), allocatable, intent(inout)    :: velocities(:,:)
    real(pr)                                :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities)/(3.0_pr * real(num_atoms,pr))

    scaling_factor = sqrt( 1._pr + dt/Berendsen_time*(ref_Temp / instant_Temp -1._pr))

    velocities = velocities*scaling_factor

end subroutine thermostat_Berendsen

END MODULE thermostatsMod