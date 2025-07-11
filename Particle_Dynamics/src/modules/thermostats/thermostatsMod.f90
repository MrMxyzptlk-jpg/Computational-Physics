MODULE thermostatsMod
    use parametersMod
    use propertiesMod
    use MD_rescaleMod
    use MD_BerendsenMod
    use MC_RndStepMod
    implicit none

    abstract interface
        subroutine thermo()
        end subroutine thermo
    end interface

    procedure(thermo), pointer     :: thermostat_chosen   => null()

END MODULE thermostatsMod