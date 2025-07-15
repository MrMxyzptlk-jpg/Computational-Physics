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

    MC_delta = min(MC_delta, 0.5_pr*sigma)
    if (MC_delta==0.5_pr*sigma) print *, "WARNING: MC_delta reached cap value sigma/2"
    print*, "MC_delta = ", MC_delta, "  Accepted ratio = ", real(MC_accepted,pr)/real(num_atoms*thermostat_steps,pr) &
        ,"  Accepted trials = ", MC_accepted
    MC_accepted = 0

end subroutine update_random_step

END MODULE MC_RndStepMod