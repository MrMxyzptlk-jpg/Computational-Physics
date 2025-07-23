! Module with all subroutines to write to files. Should be separated to different modules for easier maintenance.
MODULE writing2filesMod
    use precisionMod
    use dimensionsMod
    use subroutinesMod
    use formatsMod
    use parsingMod
    use propertiesMod, only: symbols, positions, velocities, charges, dipoles
    use FoX_dom
    implicit none

    private     size_x, size_y, size_z, symbol, dataDir
    character(len=11)       :: size_x, size_y, size_z
    character(len=5)       :: dataDir = "data/"
    character(len=3)        :: symbol

    integer(int_medium)     :: unit_positions, unit_observables, unit_structFact

contains

subroutine initialize_XYZ_data()

    write (size_x,'(E11.5)') redimensionalize(periodicity(1), "distance")
    write (size_y,'(E11.5)') redimensionalize(periodicity(2), "distance")
    write (size_z,'(E11.5)') redimensionalize(periodicity(3), "distance")

    size_x = adjustl(trim(size_x))
    size_y = adjustl(trim(size_y))
    size_z = adjustl(trim(size_z))

end subroutine initialize_XYZ_data

subroutine open_XYZ_file()
    if (save_positions)  open(newunit=unit_positions, file=dataDir//"positions.xyz", status="replace")
end subroutine open_XYZ_file

subroutine open_observables_file(reciprocal_vec)
    real (pr), intent (in)  :: reciprocal_vec(3)

    if (save_observables) then
        open(newunit=unit_observables, file=dataDir//"observables.out", status="replace")
        if (integrator == 'velocity-Verlet') then
            write(unit_observables,*) "##    t[s]    |    e_pot    |     e_kin      |   Pressure    |   Temperature"
        else
            write(unit_observables,*) "##    t[s]    |    e_pot    |     Pressure    "
        end if
    end if

    if (do_structure_factor) then
        open(newunit=unit_structFact, file=dataDir//"structure_factor.out", status="replace")
        write(unit_structFact, '(a,3(E12.5,1x))') "## Reciprocal vector: K = ", reciprocal_vec
        write(unit_structFact, '(a,2(I2,a),I2)')  "## Miller indexes: h = ", Miller_index(1), " ; k = ", Miller_index(2) &
            , " ; l = ", Miller_index(3)
        write(unit_structFact, '(a)') "## t | structure_factor(K,t)"
    end if
end subroutine open_observables_file

subroutine write_tasks(time, energies, pressures, temperatures, structure_factor)
    real (pr), intent (in)                  :: time, energies(2), pressures, temperatures, structure_factor

    if (save_observables)    call write_observables(time, energies, pressures, temperatures)
    if (do_structure_factor) call write_structure_factor(time, structure_factor)

end subroutine write_tasks

subroutine close_XYZ_file()
    if (save_positions)         close (unit_positions)
end subroutine close_XYZ_file

subroutine close_observables_file()
    if (save_observables)       close (unit_observables)
    if (do_structure_factor)    close (unit_structFact)
end subroutine close_observables_file

subroutine write_XYZfile(time)
    real (pr), intent (in)              :: time
    integer                             :: i
    character(len=11)                   :: time_tmp

    symbol = 'X'      ! Change to real element if needed

    ! Write time in string format
    write(time_tmp,'(E11.5)') redimensionalize(time, "time")
    time_tmp = adjustl(trim(time_tmp))

    ! Line 1: number of particles
    write(unit_positions, '(i6)') num_atoms

    if (integrator /= 'velocity-Verlet') then      ! Monte Carlo method does not consider velocities
        ! Line 2: extended XYZ header with box info and time
        write(unit_positions,'(A)') 'Lattice="' // &
            size_x // ' 0.0  0.0  0.0 ' // &
            size_y // ' 0.0  0.0  0.0 ' // &
            size_z // '" Properties=species:S:1:pos:R:3 Time=' // &
            time_tmp

        do i = 1, num_atoms
            if (allocated(symbols)) then
                write(unit_positions, fmt=format_XYZ) symbols(i), redimensionalize(positions(:,i), "distance")
            else
                write(unit_positions, fmt=format_XYZ) symbol, redimensionalize(positions(:,i), "distance")
            end if
        end do
    else if (integrator == 'velocity-Verlet') then
        ! Line 2: extended XYZ header with box info and time
        write(unit_positions,'(A)') 'Lattice="' // &
            size_x // ' 0.0  0.0  0.0 ' // &
            size_y // ' 0.0  0.0  0.0 ' // &
            size_z // '" Properties=species:S:1:pos:R:3:vel:R:3 Time=' // &
            time_tmp

        do i = 1, num_atoms
            if (allocated(symbols)) then
                write(unit_positions,fmt=format_XYZ) symbols(i), redimensionalize(positions(:,i), "distance") &
                    , redimensionalize(velocities(:,i), "velocity")
            else
                write(unit_positions,fmt=format_XYZ) symbol, redimensionalize(positions(:,i), "distance") &
                    , redimensionalize(velocities(:,i), "velocity")
            end if
        end do
    end if

end subroutine write_XYZfile

subroutine write_stateXML()
    character(len=9)    :: filename = 'STATE.xml'
    type(Node), pointer :: doc, root, atomsNode, atomNode, physicalNode
    character(len=128)   :: attr_string
    type(Node), pointer :: dummy
    integer             :: i

    ! Create a new XML document
    doc => createDocument(getImplementation(), "", "STATE", null())
    root => getDocumentElement(doc)

    ! Create <atoms> node with num_atoms and periodicity attributes
    physicalNode => createElementNS(doc, "", "physical")

    write(attr_string, '(I6)') num_atoms
    call setAttribute(physicalNode, "num_atoms", trim(adjustl(attr_string)))

    write(attr_string, format_state) redimensionalize(periodicity, "distance")
    call setAttribute(physicalNode, "periodicity", trim(adjustl(attr_string)))

    write(attr_string, '(a)') structure
    call setAttribute(physicalNode, "structure", trim(adjustl(attr_string)))

    dummy => appendChild(root, physicalNode)

    ! Write atoms' nodes
    atomsNode => createElementNS(doc, "", "atoms")
    dummy => appendChild(root, atomsNode)

    do i = 1, num_atoms
        atomNode => createElementNS(doc, "", "atom")

        ! Symbol (optional)
        if (allocated(symbols)) then
            call setAttribute(atomNode, "symbol", trim(symbols(i)))
        end if

        ! Position
        write(attr_string, format_state) redimensionalize(positions(:,i), "distance")
        call setAttribute(atomNode, "position", trim(adjustl(attr_string)))

        ! Velocity (optional)
        if (allocated(velocities)) then
            write(attr_string, format_state) redimensionalize(velocities(:,i), "velocity")
            call setAttribute(atomNode, "velocity", trim(adjustl(attr_string)))
        end if

        ! Charge (optional)
        if (allocated(charges)) then
            write(attr_string, *) charges(i)
            call setAttribute(atomNode, "charge", trim(adjustl(attr_string)))
        end if

        ! Dipole (optional)
        if (allocated(dipoles)) then
            write(attr_string, format_state) dipoles(1,i), dipoles(2,i), dipoles(3,i)
            call setAttribute(atomNode, "dipole", trim(adjustl(attr_string)))
        end if

        dummy => appendChild(atomsNode, atomNode)
    end do

    ! Write to file
    call serialize(doc, filename)
    print*, "STATE.xml file will not be pretty but can be fixed with 'xmllint --format STATE.xml --output STATE.xml'"

    ! Free memory
    call destroy(doc)

end subroutine write_stateXML

subroutine write_observables(time, energies, pressures, temperatures)
    real (pr), intent (in)  :: time, energies(2), pressures, temperatures

    select case(integrator)
        case ('velocity-Verlet')
            write(unit_observables, format_style0) redimensionalize(time, "time"), redimensionalize(energies, "energy") &
                , redimensionalize(pressures, "pressure"), temperatures
        case('Monte-Carlo')
            write(unit_observables, format_style0) redimensionalize(time, "time"), redimensionalize(energies(1), "energy") &
                , redimensionalize(pressures, "pressure")
        case('Brownian')
            write(unit_observables, format_style0) redimensionalize(time, "time"), redimensionalize(energies(1), "energy") &
                , redimensionalize(pressures, "pressure")
        case default
            write(unit_observables, format_style0) redimensionalize(time, "time"), redimensionalize(energies(1), "energy") &
                , redimensionalize(pressures, "pressure")
    end select

end subroutine write_observables

subroutine write_structure_factor(time, structure_factor)
    real (pr), intent (in)  :: structure_factor, time

    write(unit_structFact, format_style0) redimensionalize(time, "time"), structure_factor

end subroutine write_structure_factor

subroutine write_msd(msd)
    real (pr), intent (in)  :: msd(0:)
    real (pr)               :: time_jump
    integer                 :: unitnum, i

    time_jump = redimensionalize(dt, "time")*measuring_jump

    open(newunit=unitnum, file=dataDir//"mean_sqr_displacement.out", status="replace")
        write(unitnum, '(a)') "## Δt | ⟨Δr(t)^2⟩"
        do i = 0, size(msd) - 1
            write(unitnum, format_style0) real(i,pr)*time_jump, redimensionalize(msd(i), "distance")
        end do
    close(unitnum)

end subroutine write_msd

subroutine write_pair_corr(pair_corr)
    real (pr), intent (in)  :: pair_corr(:)
    real (pr)               :: bin_center
    integer                 :: unitnum, i

    open(newunit=unitnum, file=dataDir//"pair_correlation.out", status="replace")
        write(unitnum, '(a)') "## r | pair_correlation(r)"
        do i = 1, pair_corr_bins
            bin_center = (real(i-1,pr) + 0.5_pr)*redimensionalize(dr, "distance")
            write(unitnum, format_style0) bin_center, pair_corr(i)
        end do
    close(unitnum)

end subroutine write_pair_corr

subroutine write_output(CPU_elapsed_time, energies, pressures, temperatures, structure_factor)
    real (pr), intent (in)  :: CPU_elapsed_time, energies(:,:), pressures(:), temperatures(:), structure_factor(:)
    real (pr)               :: avg, E_total(size(energies,2))
    real (pr)               :: stddev, error
    integer                 :: unit_info

    open(newunit=unit_info, file=dataDir//"INFO.out", status="replace")
        write(unit_info,'(a24,13x,I4)')     "Number of atoms:       ", num_atoms
        write(unit_info,'(a24,14x,a)')      "Initial structure:     ", structure
        write(unit_info,'(a24,6x,E11.5)')   "Lattice constant:      ", redimensionalize(lattice_constant, "distance")
        write(unit_info,'(a24,5x,3(x,E11.5))')   "Super-cell dimensions: ", redimensionalize(periodicity, "distance")
        write(unit_info,'(a24,6x,E11.5)')   "Density:               ", redimensionalize(density, "density")
        if (interactions == 'Coulomb') write(unit_info,'(a24,6x,E11.5)')   "Gamma:                 " &
            , ((fourPi/3._pr)*redimensionalize(density, "density"))**(1._pr/3._pr)/(redimensionalize(ref_Temp, "temperature"))
        write(unit_info,'(a24,4x,a)')       "Summation:             ", summation
        write(unit_info,'(a24,4x,a)')       "Initial velocities:    ", initial_velocities
        write(unit_info,'(a24,4x,a)')       "Potential:             ", interactions
        write(unit_info,'(a24,2x,a)')       "Integrator:            ", integrator
        write(unit_info,'(a24,11x,I6)')     "Transitory steps:      ", transitory_steps
        write(unit_info,'(a24,11x,I6)')     "Run steps:             ", real_steps
        write(unit_info,'(a24,11x,I6)')     "Measuring steps:       ", measuring_steps
        write(unit_info,'(a24,6x,E11.5)')   "Simulated time:        ", real(real_steps,pr)*redimensionalize(dt, "time")

        if (save_observables) then
            select case(integrator)
            case ('velocity-Verlet')
                call get_stats(redimensionalize(pressures(1:), "pressure"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Pressure:          ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                call get_stats(redimensionalize(temperatures(1:), "temperature"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Temperature:       ", " Average = ",avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                call get_stats(redimensionalize(energies(1,1:), "energy"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Potential Energy:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                call get_stats(redimensionalize(energies(2,1:), "energy"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Kinetic Energy:    ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                E_total = energies(1,:) + energies(2,:)
                call get_stats(redimensionalize(E_total, "energy"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Total Energy:      ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error

                if (do_structure_factor) then
                    call get_stats(structure_factor, average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Structure Factor:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                end if

            case('Monte-Carlo')
                call get_stats(redimensionalize(energies(1,1:), "energy"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Potential Energy:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                call get_stats(redimensionalize(pressures(1:), "pressure"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Pressure:          ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                if (do_structure_factor) then
                    call get_stats(structure_factor, average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Structure Factor:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                end if

            case('Brownian')
                call get_stats(redimensionalize(energies(1,1:), "energy"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Potential Energy:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                call get_stats(redimensionalize(pressures(1:), "pressure"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Pressure:          ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                if (do_structure_factor) then
                    call get_stats(structure_factor, average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Structure Factor:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                end if

            case default
                call get_stats(redimensionalize(energies(1,1:), "energy"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Potential Energy:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                call get_stats(redimensionalize(pressures(1:), "pressure"), average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Pressure:          ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                if (do_structure_factor) then
                    call get_stats(structure_factor, average = avg, stddev = stddev, error = error)
                write(unit_info,format_observables)     "Structure Factor:  ", " Average = ", avg &
                    , " Standard deviation = ", stddev, " Standard error = ", error
                end if
            end select
        end if
        write(unit_info,'(a)')"-----------------------------------------------------------------------------"
        write(unit_info,'(a20,5x,F11.5, a)')"Elapsed time:          ", CPU_elapsed_time,"s"
    close(unit_info)

end subroutine write_output

END MODULE writing2filesMod