MODULE phase_initializeMod
    use initializationsMod
    use parsingMod
    use thermostatsMod
    use writing2filesMod
    implicit none

    character(len=3)    :: parser

CONTAINS

subroutine phase_initialize()

    !###################################################################################################
    !   Set default values and parse input file (parsing modules)
    !###################################################################################################

        call set_defaults()

        print*, "Parser (NML or XML): "
        read(*,*) parser
        select case(parser)
            case("NML")
                call parse_inputNML()
            case("XML")
                call parse_inputXML()
            case default
                print*, "ERROR: invalid parser" ; STOP
        end select

        call check_inputValues()

    !##################################################################################################
    !      Necessary definitions, pointers, initializations and conversion factors (all in initializationsMod module unless specified otherwise)
    !##################################################################################################

        call init_structure()
        call init_potential()
        call init_variables()
        call init_positions()
        call init_thermostat()
        call init_observables()
        call init_tasks()
        call init_integrator()
        call init_summation()
        call init_internal_constants()

        if (integrator == 'velocity-Verlet') then
            if (state == 'fromScratch') call init_velocities()
            call thermostat_rescale()     ! thermostatsMod module
        end if

        if (summation == "Ewald") then
            call init_Ewald()
            if ((integrator == "Monte-Carlo")) call init_reciprocalCharges()
        end if

        call initialize_XYZ_data()                  ! writing2fliesMod module

end subroutine phase_initialize

END MODULE phase_initializeMod