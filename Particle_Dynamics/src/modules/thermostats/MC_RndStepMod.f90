MODULE MC_RndStepMod
    use variablesMod
    use propertiesMod
    implicit none

CONTAINS

subroutine update_random_step()

    if (real(MC_accepted,pr)/real(num_atoms*thermostat_steps,pr) > MC_acceptance_rate) then
        MC_delta = MC_delta*1.05_pr
    else
        MC_delta = MC_delta*0.95_pr
    endif
    print*, MC_delta, real(MC_accepted,pr)/real(num_atoms*thermostat_steps,pr), MC_accepted, MC_acceptance_rate
    MC_accepted = 0

end subroutine update_random_step

END MODULE MC_RndStepMod