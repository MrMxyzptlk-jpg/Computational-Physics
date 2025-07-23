MODULE MC_RndStepMod
    use variablesMod
    use propertiesMod
    implicit none

CONTAINS

subroutine update_random_step() ! Flags should be added to prevent illogical MC_deltaMin and MC_deltaMax values

    if (real(MC_accepted,pr)/real(num_atoms*thermostat_steps,pr) > MC_acceptance_rate) then
        MC_delta = MC_delta*1.05_pr
    else
        MC_delta = MC_delta*0.95_pr
    endif

    MC_delta = min(MC_delta, MC_deltaMax)
    if (MC_delta == MC_deltaMax) print *, "WARNING: MC_delta reached maximum allowed value MC_deltaMax =", MC_deltaMax

    MC_delta = max(MC_delta, MC_deltaMin)
    if (MC_delta == MC_deltaMin) print *, "WARNING: MC_delta reached minimum allowed value MC_deltaMin", MC_deltaMin


    if (debugg) print*, "MC_delta = ", MC_delta, "  Accepted ratio = ", real(MC_accepted,pr)/real(num_atoms*thermostat_steps,pr) &
                ,"  Accepted trials = ", MC_accepted
    MC_accepted = 0

end subroutine update_random_step

END MODULE MC_RndStepMod