MODULE checkParsingMod
    use precisionMod
    use parametersMod
    use FoX_dom
    implicit none

CONTAINS

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
        else if ((checks_count(2) == num_atoms) .and. (initial_velocities == 'fromFile')) then
            checks(2) = .True.
        end if
    end if

    if (interactions == 'Coulomb') then
        ! Charges check
        if ((checks_count(3) /= 0) .and. (checks_count(3) /= num_atoms)) then
            print*, "Mismatch between 'charges' and number of atoms. Stopping..."
            STOP
        else if (checks_count(3) == 0) then
            print*, "No charges specified for 'Coulomb' interactions. Stopping..."; STOP
        else if (checks_count(3) == num_atoms) then
            checks(3) = .True.
        end if

        ! Dipoles check
        !if ((checks_count(4) /= 0) .and. (checks_count(4) /= num_atoms)) then
        !    print*, "Mismatch between 'dipoles' and number of atoms. Stopping..."
        !    STOP
        !else if (checks_count(4) == 0) then
        !    print*, "No dipoles specified for 'Coulomb' interactions. Stopping..."; STOP
        !else if (checks_count(4) == num_atoms) then
        !    checks(4) = .True.
        !end if
    end if

end subroutine check_stateProperties

subroutine check_inputValues()
    call check_parsed_physical()
    call check_parsed_calculation()
    call check_parsed_tasks()
    call check_parsed_approximation()
    call check_parsed_thermostat()
    if (do_mean_sqr_displacement) call check_parsed_MSD()
    if (do_pair_correlation) call check_parsed_pair_correlation()
end subroutine check_inputValues

!############ Subroutines to check values of each parsed section ###############

subroutine check_parsed_physical()
    if ((ensemble /= 'NVE') .and. (ensemble /= 'NVT')) print*, "ERROR: unavailable ensemble"; STOP
    if ((structure /= 'random') .and. (structure /= 'FCC') .and. (structure /= 'BCC')) print*, "ERROR: unavailable structure"; STOP
    if (lattice_constant <= 0._pr) print*, "ERROR: lattice_constant <= 0"; STOP
    if (cell_dim(1) <= 0._pr) print*, "ERROR: 1st value of cell_dim <= 0"; STOP
    if (cell_dim(2) <= 0._pr) print*, "ERROR: 2nd value of cell_dim <= 0"; STOP
    if (cell_dim(3) <= 0._pr) print*, "ERROR: 3rd value of cell_dim <= 0"; STOP
    if (num_atoms <= 0) print*, "ERROR: num_atoms <= 0"; STOP
    if (ref_Temp <= 0._pr) print*, "ERROR: ref_Temp <= 0"; STOP
    if (density <= 0._pr) print*, "ERROR: density <= 0"; STOP
    if (viscosity <= 0._pr) print*, "ERROR: viscosity <= 0"; STOP
    if (reduced_viscosity   <= 0._pr) print*, "ERROR: reduced_viscosity <= 0"; STOP
    if (mass <= 0._pr) print*, "ERROR: mass <= 0"; STOP
end subroutine check_parsed_physical

subroutine check_parsed_calculation()
    if ((state /= 'fromScratch') .and. (state /= 'fromFile')) print*, "ERROR: invalid state option"; STOP
    if (real_steps <= 0) print*, "ERROR: real_steps <= 0"; STOP
    if (transitory_steps <= 0) print*, "ERROR: transitory_steps <= 0"; STOP
    if (thermostat_steps <= 0) print*, "ERROR: thermostat_steps <= 0"; STOP
    if (dt <= 0) print*, "ERROR: dt <= 0"; STOP
    if (radius_cutoff <= 0) print*, "ERROR: radius_cutoff <= 0"; STOP
    if (measuring_jump <= 0) print*, "ERROR: measuring_jump <= 0"; STOP
    if ((initial_velocities /= 'random') .and. (initial_velocities /= 'Maxwell') .and. (initial_velocities /= 'fromfile')) print*, &
        "ERROR: invalid initial_velocities"; STOP
    if ((summation /= 'all-vs-all') .and. (summation /= 'linked-lists') .and. (summation /= 'Ewald')) print*, &
        "ERROR: unavailable summation"; STOP
    if (summation == 'linked-lists') then   ! Another check is done once the periodicity is calculated
        if (dim_linkCell(1) <= 0) print*, "ERROR: 1st value of dim_linkCell <= 0"; STOP
        if (dim_linkCell(2) <= 0) print*, "ERROR: 2nd value of dim_linkCell <= 0"; STOP
        if (dim_linkCell(3) <= 0) print*, "ERROR: 3rd value of dim_linkCell <= 0"; STOP
    end if
end subroutine check_parsed_calculation

subroutine check_parsed_tasks()
    if ((save_transitory .neqv. .True.) .and. (save_transitory .neqv. .False.)) print*, "ERROR: invalid save_transitory"; STOP
    if ((save_observables .neqv. .True.) .and. (save_observables .neqv. .False.)) print*, "ERROR: invalid save_observables"; STOP
    if ((save_positions .neqv. .True.) .and. (save_positions .neqv. .False.)) print*, "ERROR: invalid save_positions"; STOP
    if ((save_state .neqv. .True.) .and. (save_state .neqv. .False.)) print*, "ERROR: invalid save_state"; STOP
    if ((do_pair_correlation .neqv. .True.) .and. (do_pair_correlation .neqv. .False.)) &
        print*, "ERROR: invalid do_pair_correlation"; STOP
    if ((do_mean_sqr_displacement .neqv. .True.) .and. (do_mean_sqr_displacement .neqv. .False.)) &
        print*, "ERROR: invalid do_mean_sqr_displacement"; STOP
    if ((do_structure_factor .neqv. .True.) .and. (do_structure_factor .neqv. .False.)) &
        print*, "ERROR: invalid do_structure_factor"; STOP
end subroutine check_parsed_tasks

subroutine check_parsed_approximation()
    if ((integrator /= 'velocity-Verlet') .and. (integrator /= 'Brownian') .and. (integrator /= 'Monte-Carlo')) then
        print*, "ERROR: unavailable integrator"
        STOP
    end if

    if (sigma <= 0._pr)     print*, "ERROR: sigma <= 0._pr"; STOP
    if (epsilon <= 0._pr)   print*, "ERROR: epsilon <= 0._pr"; STOP

    if (integrator == 'Monte-Carlo') then
        if (MC_adjust_step <= 0) print*, "ERROR: MC_adjust_step <= 0"; STOP
        if (MC_delta <= 0._pr)   print*, "ERROR: MC_delta <= 0._pr"; STOP
    end if

    if  ((integrator == 'velocity-Verlet') .and. (interactions == 'coulomb') .and. (summation /= 'Ewald')) then
        print*, "WARNING: unavailable summation for Coulomb interactions. Switching to 'Ewald' summation"
        summation = 'Ewald'
    end if

    if ((integrator == 'velocity-Verlet') .and. (summation == 'Ewald')) then
        if (kgrid(1) <= 0) print*, "ERROR: 1st value of kgrid <= 0"; STOP
        if (kgrid(2) <= 0) print*, "ERROR: 2nd value of kgrid <= 0"; STOP
        if (kgrid(3) <= 0) print*, "ERROR: 3rd value of kgrid <= 0"; STOP
    end if
end subroutine check_parsed_approximation

subroutine check_parsed_thermostat()
    if ((thermostat_type /= 'rescale') .and. (thermostat_type /= 'Berendsen')) then
        print*, "ERROR: unavailable thermostat", thermostat_type
        STOP
    end if

    if ((thermostat_type == 'Berendsen').and.(Berendsen_time <= 1e-10)) then
        print*, "ERROR: Berendsen_time <= 1e-10"
        STOP
    end if
end subroutine check_parsed_thermostat

subroutine check_parsed_MSD()
    if (max_correlation <= 0) print*, "ERROR: max_correlation <= 0"; STOP
    if (max_correlation <= real_steps/measuring_jump) print*, "WARNING: max_correlation <= real_steps/measuring_jump"
end subroutine check_parsed_MSD

subroutine check_parsed_pair_correlation()
    ! pair_corr_cutoff must be checked as soon as the periodicity is calculated
    if (pair_corr_bins <= 10) print*, "ERROR: pair_corr_bins <= 10"; STOP
end subroutine check_parsed_pair_correlation

END MODULE checkParsingMod