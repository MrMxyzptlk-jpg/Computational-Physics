MODULE parsingMod
    use precisionMod
    use subroutinesMod
    use FoX_dom

    implicit none

    private inputNode
    type(Node), pointer     :: inputNode

    character (len=15)  :: type
    character (len=6)   :: ensemble
    character (len=12)  :: summation, thermostat_type, initial_velocities
    logical             :: save_positions, do_structure_factor, do_mean_sqr_displacement
    integer(int_large)  :: transitory_steps, thermostat_steps, dim_linkCell(3), Miller_index(3)

    interface get_parsed_value
        module procedure get_parsed_string
        module procedure get_parsed_logical
        module procedure get_parsed_integerScalar
        module procedure get_parsed_integerVector
        module procedure get_parsed_realScalar
        module procedure get_parsed_realVector
    end interface


    ! Namelist blocks
    namelist /physical/ structure, lattice_constant, density, reduced_viscosity, viscosity, ref_Temp, num_atoms, mass, cell_dim &
        , ensemble
    namelist /calculation/ real_steps, transitory_steps, thermostat_steps, dt, radius_cutoff, summation, dim_linkCell &
        , measuring_jump, initial_velocities
    namelist /tasks/ save_transitory, save_positions, save_observables, do_pair_correlation, do_mean_sqr_displacement &
        , do_structure_factor, Miller_index
    namelist /approximation/ integrator, type, sigma, epsilon, kgrid, MC_adjust_step, MC_delta
    namelist /thermostat/ thermostat_type, Berendsen_time
    namelist /MSD/ max_correlation
    namelist /pair_correlation/pair_corr_cutoff, pair_corr_bins

    CONTAINS

subroutine set_defaults()
        ! Physical problems' characteristics
        ensemble            = "NVE"
        structure           = "FCC"
        lattice_constant    = 1._pr
        cell_dim            = (/1_int_small,1_int_small,1_int_small/)
        num_atoms           = 0
        ref_Temp            = 1.0_pr
        density             = 0._pr
        viscosity           = 0._pr
        reduced_viscosity   = 0._pr
        mass                = 1._pr

        !Calculation settings
        real_steps       = 1000
        transitory_steps = 1000
        thermostat_steps = 50
        dt               = 0.005_pr
        radius_cutoff    = 2.5_pr
        summation        = 'all-vs-all'
        dim_linkCell     = (/3,3,3/)
        measuring_jump   = 50
        initial_velocities = 'random'

        ! Tasks
        save_transitory  = .False.
        save_observables = .False.
        save_positions   = .False.
        do_pair_correlation = .False.
        do_mean_sqr_displacement = .False.
        do_structure_factor = .False.
        Miller_index        = (/0,-1,0/)

        ! Potential parameters
        integrator  = 'velocity-Verlet'
        sigma   = 1._pr
        epsilon = 1._pr
        kgrid   = (/5, 5, 5/)
        MC_adjust_step  = 50
        MC_delta        = 0.01

        ! Thermostat parameters
        thermostat_type      = 'rescale'
        Berendsen_time  = 0.01_pr

        ! MSD
        max_correlation  = 100

        ! Pair correlation parameters
        pair_corr_cutoff = 4.0_pr
        pair_corr_bins   = 100

end subroutine set_defaults

subroutine parse_inputNML()
    integer                 :: unit_input
    ! Read from input file
    open(newunit=unit_input, file="input.nml", status="old", action="read")
        read(unit_input, nml=physical)
        read(unit_input, nml=calculation)
        read(unit_input, nml=tasks)
        read(unit_input, nml=approximation)

        if(ensemble=="NVT") then
            read(unit_input, nml=thermostat)
        end if

        if(do_mean_sqr_displacement) then
            read(unit_input, nml=MSD)
        end if

        if(do_pair_correlation) then
            read(unit_input, nml=pair_correlation)
        end if

    close(unit_input)

end subroutine parse_inputNML

subroutine parse_inputXML()
    character(len=9)        :: inputFile = 'input.xml'
    type(Node), pointer     :: inputDoc
    type(NodeList), pointer :: list

    call check_fileXML(inputFile, inputDoc)

    !####### PHYSICAL node #######
    list => getElementsByTagName(inputDoc, "physical")
    inputNode => item(list, 0)

    call get_parsed_value("ensemble", ensemble)
    call get_parsed_value("structure", structure)
    call get_parsed_value("lattice_constant", lattice_constant)
    call get_parsed_value("cell_dim", cell_dim)
    call get_parsed_value("num_atoms", num_atoms)
    call get_parsed_value("ref_Temp", ref_Temp)
    call get_parsed_value("density", density)
    call get_parsed_value("reduced_viscosity", reduced_viscosity)
    call get_parsed_value("mass", mass)


    !####### CALCULATION node #######
    list => getElementsByTagName(inputDoc, "calculation")
    inputNode => item(list, 0)

    call get_parsed_value("dt", dt)
    call get_parsed_value("real_steps", real_steps)
    call get_parsed_value("transitory_steps", transitory_steps)
    call get_parsed_value("thermostat_steps", thermostat_steps)
    call get_parsed_value("radius_cutoff", radius_cutoff)
    call get_parsed_value("summation", summation)
    call get_parsed_value("dim_linkCell", dim_linkCell)
    call get_parsed_value("measuring_jump", measuring_jump)
    call get_parsed_value("initial_velocities", initial_velocities)
    if (initial_velocities == "Maxwell") then
        print*, "Maxwell not debugged for XML parser. Changing to random"
        initial_velocities = "random"
    end if

    !####### TASKS node #######
    list => getElementsByTagName(inputDoc, "tasks")
    inputNode => item(list, 0)

    call get_parsed_value("save_transitory", save_transitory)
    call get_parsed_value("save_observables", save_observables)
    call get_parsed_value("save_positions", save_positions)
    call get_parsed_value("do_pair_correlation", do_pair_correlation)
    call get_parsed_value("do_structure_factor", do_structure_factor)
    call get_parsed_value("do_mean_sqr_displacement", do_mean_sqr_displacement)
    call get_parsed_value("Miller_index", Miller_index)

    !####### APPROXIMATION node #######
    list => getElementsByTagName(inputDoc, "approximation")
    inputNode => item(list, 0)

    call get_parsed_value("integrator", integrator)
    call get_parsed_value("type", type)
    call get_parsed_value("sigma", sigma)
    call get_parsed_value("epsilon", epsilon)
    call get_parsed_value("kgrid", kgrid)
    call get_parsed_value("MC_adjust_step", MC_adjust_step)
    call get_parsed_value("MC_delta", MC_delta)

    !####### THERMOSTAT node #######
    if(ensemble=="NVT") then
        list => getElementsByTagName(inputDoc, "thermostat")
        inputNode => item(list, 0)

        call get_parsed_value("thermostat_type", thermostat_type)
        call get_parsed_value("Berendsen_time", Berendsen_time)
    end if

    !####### MSD node #######
    if(do_mean_sqr_displacement) then
        list => getElementsByTagName(inputDoc, "MSD")
        inputNode => item(list, 0)

        call get_parsed_value("max_correlation", max_correlation)
    end if

    !####### pair_correlation node #######
    if(do_pair_correlation) then
        list => getElementsByTagName(inputDoc, "pair_correlation")
        inputNode => item(list, 0)

        call get_parsed_value("pair_corr_cutoff", pair_corr_cutoff)
        call get_parsed_value("pair_corr_bins", pair_corr_bins)
    end if

    ! Clean up
    call destroy(inputDoc)

end subroutine parse_inputXML

subroutine check_fileXML(filename, fileDoc)
    character(len=*), intent(in)        :: filename
    type(Node), pointer, intent(out)    :: fileDoc
    type(DOMException)                  :: ex
    integer                             :: ios

    fileDoc => parseFile(filename, iostat=ios, ex=ex)
    if (ios /= 0) then
        print*, "Error reading file. iostat was ", ios
        stop
    else if (inException(ex)) then
        print*,"DOM Parse error ", getExceptionCode(ex)
        stop
    else
        print*, "Read input.xml file."
    endif

end subroutine check_fileXML

subroutine get_parsed_string(name, value)
    character(len=*), intent(in)    :: name
    character(len=*)                :: value

    if (hasAttribute(inputNode, name)) then
        value = getAttribute(inputNode, name)
        print*, name, value
    end if

end subroutine get_parsed_string

subroutine get_parsed_logical(name, value)
    character(len=*), intent(in)    :: name
    logical                         :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
        print*, name, value
    end if

end subroutine get_parsed_logical

subroutine get_parsed_integerScalar(name, value)
    character(len=*), intent(in)    :: name
    integer                         :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
        print*, name, value
    end if

end subroutine get_parsed_integerScalar

subroutine get_parsed_integerVector(name, value)
    character(len=*), intent(in)    :: name
    integer                         :: value(3)
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value(1), value(2), value(3)
        print*, name, value
    end if

end subroutine get_parsed_integerVector

subroutine get_parsed_realScalar(name, value)
    character(len=*), intent(in)    :: name
    real(pr)                        :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
        print*, name, value
    end if

end subroutine get_parsed_realScalar

subroutine get_parsed_realVector(name, value)
    character(len=*), intent(in)    :: name
    real(pr)                        :: value(3)
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value(1), value(2), value(3)
        print*, name, value
    end if

end subroutine get_parsed_realVector


!##################################################################################################
!     Not used / Not implemented
!##################################################################################################

!!!subroutine parse_velocities(num_atoms, positions)
!!!    integer                 :: unit_input
!!!    ! Read from input file
!!!    open(newunit=unit_input, file="positions.in", status="old", action="read")
!!!
!!!    close(unit_input)
!!!
!!!end subroutine parse_velocities
!!!
!!!subroutine parse_positions(num_atoms, positions)
!!!    integer, intent(inout)  :: num_atoms
!!!    real(pr), intent(inout)   :: position(:)
!!!    character(5)            :: symbols
!!!    integer                 :: unit_input, parsed_num_atoms
!!!    ! Read from input file
!!!    open(newunit=unit_input, file="STATE.xyz", status="old", action="read")
!!!        read(*,'(I10)')   parsed_num_atoms
!!!        if (parsed_num_atoms /= num_atoms) then
!!!            print*, "Input num_atoms does not coincide with the one parsed from STATE.xyz"
!!!            print*, "Using parsed num_atoms"
!!!        end if
!!!    close(unit_input)
!!!
!!!end subroutine parse_positions

END MODULE parsingMod