MODULE parsingMod
    use precisionMod
    use subroutinesMod
    use propertiesMod
    use parametersMod
    use checkParsingMod
    use FoX_dom

    implicit none

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
        , measuring_jump, initial_velocities, state
    namelist /tasks/ save_transitory, save_positions, save_observables, do_pair_correlation, do_mean_sqr_displacement &
        , do_structure_factor, Miller_index, save_state
    namelist /approximation/ integrator, interactions, sigma, delta, kgrid, MC_adjust_step, MC_delta
    namelist /thermostat/ thermostat_type, Berendsen_time
    namelist /MSD/ max_correlation
    namelist /pair_correlation/pair_corr_cutoff, pair_corr_bins

    CONTAINS

subroutine check_statePhysical(inputNode, parsed_periodicity)   ! Should be in checkParsingMod but cannot because of the use of get_parsed_value()
    type(Node), pointer, intent(in) :: inputNode
    real(pr), intent(out)           :: parsed_periodicity(3)
    integer(int_large)              :: parsed_num_atoms
    character (len=6)               :: parsed_structure

    ! Check necessary attributes exist
    if (.not.(hasAttribute(inputNode, "num_atoms"))) STOP "ERROR: Missing num_atoms in STATE.xml"
    if (.not.(hasAttribute(inputNode, "periodicity"))) STOP "ERROR: Missing periodicity in STATE.xml"
    if (.not.(hasAttribute(inputNode, "structure"))) STOP "ERROR: Missing structure in STATE.xml"

    ! Check consistent number of atoms specified
    call get_parsed_value(inputNode, "num_atoms", parsed_num_atoms)
    if (num_atoms /= parsed_num_atoms) STOP "ERROR: Mismatch between STATE.xml num_atoms and input.xml num_atoms"

    ! Check consistent structure specification
    call get_parsed_value(inputNode, "structure", parsed_structure)
    if (structure /= parsed_structure) STOP "ERROR: Mismatch between STATE.xml structure and input.xml structure"

    ! Get the other relevant variables
    call get_parsed_value(inputNode, "periodicity", parsed_periodicity)

end subroutine check_statePhysical

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
        state            = 'fromScratch'
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
        save_state       = .False.
        do_pair_correlation = .False.
        do_mean_sqr_displacement = .False.
        do_structure_factor = .False.
        Miller_index        = (/0,-1,0/)

        ! Potential parameters
        integrator  = 'velocity-Verlet'
        sigma   = 1._pr
        delta = 1._pr
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
    type(Node), pointer     :: inputNode
    type(NodeList), pointer :: list

    call check_fileXML(inputFile, inputDoc)

    !####### PHYSICAL node #######
    list => getElementsByTagName(inputDoc, "physical")
    inputNode => item(list, 0)

    call get_parsed_value(inputNode, "ensemble", ensemble)
    call get_parsed_value(inputNode, "structure", structure)
    call get_parsed_value(inputNode, "lattice_constant", lattice_constant)
    call get_parsed_value(inputNode, "cell_dim", cell_dim)
    call get_parsed_value(inputNode, "num_atoms", num_atoms)
    call get_parsed_value(inputNode, "ref_Temp", ref_Temp)
    call get_parsed_value(inputNode, "density", density)
    call get_parsed_value(inputNode, "reduced_viscosity", reduced_viscosity)
    call get_parsed_value(inputNode, "mass", mass)


    !####### CALCULATION node #######
    list => getElementsByTagName(inputDoc, "calculation")
    inputNode => item(list, 0)

    call get_parsed_value(inputNode, "dt", dt)
    call get_parsed_value(inputNode, "real_steps", real_steps)
    call get_parsed_value(inputNode, "transitory_steps", transitory_steps)
    call get_parsed_value(inputNode, "thermostat_steps", thermostat_steps)
    call get_parsed_value(inputNode, "radius_cutoff", radius_cutoff)
    call get_parsed_value(inputNode, "summation", summation)
    call get_parsed_value(inputNode, "dim_linkCell", dim_linkCell)
    call get_parsed_value(inputNode, "measuring_jump", measuring_jump)
    call get_parsed_value(inputNode, "initial_velocities", initial_velocities)

    !####### TASKS node #######
    list => getElementsByTagName(inputDoc, "tasks")
    inputNode => item(list, 0)

    call get_parsed_value(inputNode, "save_transitory", save_transitory)
    call get_parsed_value(inputNode, "save_observables", save_observables)
    call get_parsed_value(inputNode, "save_positions", save_positions)
    call get_parsed_value(inputNode, "save_state", save_state)
    call get_parsed_value(inputNode, "do_pair_correlation", do_pair_correlation)
    call get_parsed_value(inputNode, "do_structure_factor", do_structure_factor)
    call get_parsed_value(inputNode, "do_mean_sqr_displacement", do_mean_sqr_displacement)
    call get_parsed_value(inputNode, "Miller_index", Miller_index)

    !####### APPROXIMATION node #######
    list => getElementsByTagName(inputDoc, "approximation")
    inputNode => item(list, 0)

    call get_parsed_value(inputNode, "integrator", integrator)
    call get_parsed_value(inputNode, "interactions", interactions)
    call get_parsed_value(inputNode, "sigma", sigma)
    call get_parsed_value(inputNode, "delta", delta)
    call get_parsed_value(inputNode, "kgrid", kgrid)
    call get_parsed_value(inputNode, "MC_adjust_step", MC_adjust_step)
    call get_parsed_value(inputNode, "MC_delta", MC_delta)

    !####### THERMOSTAT node #######
    if(ensemble=="NVT") then
        list => getElementsByTagName(inputDoc, "thermostat")
        inputNode => item(list, 0)

        call get_parsed_value(inputNode, "thermostat_type", thermostat_type)
        call get_parsed_value(inputNode, "Berendsen_time", Berendsen_time)
    end if

    !####### MSD node #######
    if(do_mean_sqr_displacement) then
        list => getElementsByTagName(inputDoc, "MSD")
        inputNode => item(list, 0)

        call get_parsed_value(inputNode, "max_correlation", max_correlation)
    end if

    !####### pair_correlation node #######
    if(do_pair_correlation) then
        list => getElementsByTagName(inputDoc, "pair_correlation")
        inputNode => item(list, 0)

        call get_parsed_value(inputNode, "pair_corr_cutoff", pair_corr_cutoff)
        call get_parsed_value(inputNode, "pair_corr_bins", pair_corr_bins)
    end if

    ! Clean up
    call destroy(inputDoc)

end subroutine parse_inputXML

subroutine parse_stateXML(parsed_periodicity) ! Gets normalized positions, as well as other data if present.
    real(pr), intent(out)   :: parsed_periodicity(3)
    character(len=9)        :: inputFile = 'STATE.xml'
    type(Node), pointer     :: inputDoc, child
    type(NodeList), pointer :: list, children
    type(Node), pointer     :: inputNode
    integer                 :: i
    logical                 :: checks(4)

    call check_fileXML(inputFile, inputDoc)

    list => getElementsByTagName(inputDoc, "physical")
    inputNode => item(list, 0)

    call check_statePhysical(inputNode, parsed_periodicity)

    ! Check consistent properties specified and allocate properties
    list => getElementsByTagName(inputDoc, "atoms")
    inputNode => item(list, 0)
    children => getElementsByTagName(inputNode, "atom")
    call check_stateProperties(checks, children)

    allocate(positions(3,num_atoms))
    if (checks(1)) allocate(symbols(num_atoms))
    if (checks(2)) allocate(velocities(3,num_atoms))
    if (checks(3)) allocate(charges(num_atoms))
    if (checks(4)) allocate(dipoles(3,num_atoms))

    ! Get properties' values
    do i = 1, num_atoms
        child => item(children, i-1)
        call get_parsed_value(child, "position", positions(:,i)) ! Note this is not checked because it defines the number of atoms
        if (checks(1)) call get_parsed_value(child, "symbol", symbols(i))
        if (checks(2)) call get_parsed_value(child, "velocity", velocities(:,i))
        if (checks(3)) call get_parsed_value(child, "charge", charges(i))
        if (checks(4)) call get_parsed_value(child, "dipole", dipoles(:,i))
    end do

    ! Clean up
    call destroy(inputDoc)

end subroutine parse_stateXML

!############ Subroutines to inquire about and get the value of each kind of data ###############

subroutine get_parsed_string(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    character(len=*)                :: value

    if (hasAttribute(inputNode, name)) then
        value = getAttribute(inputNode, name)
    end if

end subroutine get_parsed_string

subroutine get_parsed_logical(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    logical                         :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
    end if

end subroutine get_parsed_logical

subroutine get_parsed_integerScalar(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    integer                         :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
    end if

end subroutine get_parsed_integerScalar

subroutine get_parsed_integerVector(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    integer                         :: value(3)
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value(1), value(2), value(3)
    end if

end subroutine get_parsed_integerVector

subroutine get_parsed_realScalar(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    real(pr)                        :: value
    character(len=64)               :: attr_string

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value
    end if

end subroutine get_parsed_realScalar

subroutine get_parsed_realVector(inputNode, name, value)
    character(len=*), intent(in)    :: name
    type(Node), intent(in), pointer :: inputNode
    real(pr)                        :: value(3)
    character(len=88)               :: attr_string      ! Long enough for the total format_state specifier

    if (hasAttribute(inputNode, name)) then
        attr_string = getAttribute(inputNode, name)
        read(attr_string, *) value(1), value(2), value(3)
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