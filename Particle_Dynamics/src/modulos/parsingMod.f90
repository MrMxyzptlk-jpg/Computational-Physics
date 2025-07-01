MODULE parsingMod
    use precisionMod
    use subroutinesMod
    use propertiesMod
    use FoX_dom

    implicit none

    character (len=15)  :: type, state
    character (len=6)   :: ensemble
    character (len=12)  :: summation, thermostat_type, initial_velocities
    logical             :: save_positions, save_state, do_structure_factor, do_mean_sqr_displacement
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
        , measuring_jump, initial_velocities, state
    namelist /tasks/ save_transitory, save_positions, save_observables, do_pair_correlation, do_mean_sqr_displacement &
        , do_structure_factor, Miller_index, save_state
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
    if (initial_velocities == "Maxwell") then
        print*, "Maxwell not debugged for XML parser. Changing to random"
        initial_velocities = "random"
    end if

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
    call get_parsed_value(inputNode, "type", type)
    call get_parsed_value(inputNode, "sigma", sigma)
    call get_parsed_value(inputNode, "epsilon", epsilon)
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

subroutine parse_stateXML() ! Gets normalized positions, as well as other data if present.
    character(len=9)        :: inputFile = 'STATE.xml'
    type(Node), pointer     :: inputDoc, child
    type(NodeList), pointer :: list, children
    type(Node), pointer     :: inputNode
    real(pr)                :: parsed_periodicity(3)
    integer                 :: i
    logical                 :: checks(4)

    call check_fileXML(inputFile, inputDoc)

    list => getElementsByTagName(inputDoc, "atoms")
    inputNode => item(list, 0)
    children => getChildNodes(inputNode)

    ! Check consistent number of atoms specified
    call get_parsed_value(inputNode, "num_atoms", num_atoms)
    if (getLength(children) /= num_atoms) print*, "Warning: Mismatch between num_atoms and the number of atoms in STATE.xml file."
    num_atoms = getLength(children)

    ! Check consistent properties specified and allocate properties
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

    ! Normalize positions to rescale later
    call get_parsed_value(inputNode, "periodicity", parsed_periodicity)
    print*, parsed_periodicity
    do i = 1, num_atoms
        positions(:,i) = positions(:,i) / parsed_periodicity
    end do

    ! Clean up
    call destroy(inputDoc)

end subroutine parse_stateXML

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

subroutine check_stateProperties(checks, children)
    logical                 :: checks(4)
    integer                 :: checks_count(4)
    type(NodeList), pointer :: children
    type(Node), pointer     :: child
    integer                 :: i

    checks = .False.
    checks_count = 0

    do i = 1, num_atoms
        child => item(children, i-1)
        if (hasAttribute(child, "symbol")) then
            checks_count(1) = checks_count(1) + 1
        end if
        if (hasAttribute(child, "velocity")) then
            checks_count(2) = checks_count(2) + 1
        end if
        if (hasAttribute(child, "charge")) then
            checks_count(3) = checks_count(3) + 1
        end if
        if (hasAttribute(child, "dipole")) then
            checks_count(4) = checks_count(4) + 1
        end if
    end do

    ! Symbols check
    if (checks_count(1) /= num_atoms) then
        print*, "Mismatch between 'symbols' and number of atoms. Using default symbol = X"
        checks(1) = .False.
    else if (checks_count(1) == 0) then
        checks(1) = .False.
    else if (checks_count(1) == num_atoms) then
        checks(1) = .True.
    end if

    ! Velocities check
    if (integrator == 'velocity-Verlet') then
        if ((checks_count(2) /= 0) .and. (checks_count(2) /= num_atoms)) then
            print*, "Mismatch between 'velocities' and number of atoms. Stopping..."
            STOP
        else if (checks_count(2) == 0) then
            print*, "No velocities specified for 'velocity-Verlet'. Stopping..."; STOP
        else if (checks_count(2) == num_atoms) then
            checks(2) = .True.
        end if
    end if

    if (type == 'Coulomb') then
        ! Charges check
        if ((checks_count(3) /= 0) .and. (checks_count(3) /= num_atoms)) then
            print*, "Mismatch between 'charges' and number of atoms. Stopping..."
            STOP
        else if (checks_count(3) == 0) then
            print*, "No charges specified for 'Coulomb' type interaction. Stopping..."; STOP
        else if (checks_count(3) == num_atoms) then
            checks(3) = .True.
        end if

        ! Dipoles check
        !if ((checks_count(4) /= 0) .and. (checks_count(4) /= num_atoms)) then
        !    print*, "Mismatch between 'dipoles' and number of atoms. Stopping..."
        !    STOP
        !else if (checks_count(4) == 0) then
        !    print*, "No dipoles specified for 'Coulomb' type interaction. Stopping..."; STOP
        !else if (checks_count(4) == num_atoms) then
        !    checks(4) = .True.
        !end if
    end if

end subroutine check_stateProperties

!############ Subroutines to verify and get the value of each kind of data ###############

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
    character(len=64)               :: attr_string

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