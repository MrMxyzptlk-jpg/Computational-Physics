MODULE thermostatsMod
    use parametersMod
    use propertiesMod
    implicit none
CONTAINS

subroutine thermostat_rescale()
    real(pr)    :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities)/(3.0_pr * real(num_atoms,pr))
    scaling_factor = sqrt( ref_Temp / instant_Temp )

    velocities = velocities*scaling_factor

end subroutine thermostat_rescale

subroutine thermostat_Berendsen()
    real(pr)    :: instant_Temp, scaling_factor

    instant_Temp = sum(velocities*velocities)/(3.0_pr * real(num_atoms,pr))

    scaling_factor = sqrt( 1._pr + dt/Berendsen_time*(ref_Temp / instant_Temp -1._pr))

    velocities = velocities*scaling_factor

end subroutine thermostat_Berendsen

subroutine update_random_step()
    integer, save   :: N_accepted = 0

    if (real(N_accepted,pr)/real(MC_adjust_step,pr) > 0.5_pr) then
        MC_delta = MC_delta*1.05_pr
    else
        MC_delta = MC_delta*0.95_pr
    endif
    N_accepted = 0

end subroutine update_random_step

END MODULE thermostatsMod