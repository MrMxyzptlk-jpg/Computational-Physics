MODULE MC_RndStepMod
    use variablesMod
    use propertiesMod
    implicit none

CONTAINS

subroutine update_random_step()

    if (real(MC_accepted,pr)/real(MC_adjust_step,pr) > 0.5_pr) then
        MC_delta = MC_delta*1.05_pr
    else
        MC_delta = MC_delta*0.95_pr
    endif
    MC_accepted = 0

end subroutine update_random_step

END MODULE MC_RndStepMod