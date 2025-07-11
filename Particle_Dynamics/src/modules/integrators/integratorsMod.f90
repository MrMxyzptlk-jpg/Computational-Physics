MODULE integratorsMod
    use precisionMod
    use subroutinesMod
    use parsingMod
    use updatePositionsMod
    use propertiesMod
    use forcesMod, only : get_forces, create_links
    use velVerlet_stepMod
    use MC_stepMod
    use BD_stepMod
    implicit none

    abstract interface
        subroutine integrate (step)
            use precisionMod
            integer(int_large), intent(in) :: step
        end subroutine integrate
    end interface
    procedure(integrate), pointer   :: integrator_step  => null()

END MODULE integratorsMod