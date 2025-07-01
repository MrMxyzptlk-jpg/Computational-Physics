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
    if ((ensemble /= 'NVE') .and. (ensemble /= 'NVT')) then
        print*, "ERROR: unavailable ensemble"
        STOP
    end if
    if ((structure /= 'random') .and. (structure /= 'FCC') .and. (structure /= 'BCC')) STOP "ERROR: unavailable structure"
    if (lattice_constant <= 0._pr) STOP "ERROR: lattice_constant <= 0"
    if (cell_dim(1) <= 0._pr) STOP "ERROR: 1st value of cell_dim <= 0"
    if (cell_dim(2) <= 0._pr) STOP "ERROR: 2nd value of cell_dim <= 0"
    if (cell_dim(3) <= 0._pr) STOP "ERROR: 3rd value of cell_dim <= 0"
    if (num_atoms <= 0) STOP "ERROR: num_atoms <= 0"
    if (ref_Temp <= 0._pr) STOP "ERROR: ref_Temp <= 0"
    if (density <= 0._pr) STOP "ERROR: density <= 0"
    if (mass <= 0._pr) STOP "ERROR: mass <= 0"
    ! These two are intentionally checked during initializations
    ! if (viscosity <= 0._pr) STOP "ERROR: viscosity <= 0"
    ! if (reduced_viscosity   <= 0._pr) STOP "ERROR: reduced_viscosity <= 0"
end subroutine check_parsed_physical

subroutine check_parsed_calculation()
    if ((state /= 'fromScratch') .and. (state /= 'fromFile')) STOP "ERROR: invalid state option"
    if (real_steps <= 0) STOP "ERROR: real_steps <= 0"
    if (transitory_steps <= 0) STOP "ERROR: transitory_steps <= 0"
    if (thermostat_steps <= 0) STOP "ERROR: thermostat_steps <= 0"
    if (dt <= 0) STOP "ERROR: dt <= 0"
    if (radius_cutoff <= 0) STOP "ERROR: radius_cutoff <= 0"
    if (measuring_jump <= 0) STOP "ERROR: measuring_jump <= 0"
    if ((initial_velocities /= 'random') .and. (initial_velocities /= 'Maxwell') .and. (initial_velocities /= 'fromfile')) &
        STOP "ERROR: invalid initial_velocities"
    if ((summation /= 'all-vs-all') .and. (summation /= 'linked-lists') .and. (summation /= 'Ewald')) &
        STOP "ERROR: unavailable summation"
    if (summation == 'linked-lists') then   ! Another check is done once the periodicity is calculated
        if (dim_linkCell(1) <= 0) STOP "ERROR: 1st value of dim_linkCell <= 0"
        if (dim_linkCell(2) <= 0) STOP "ERROR: 2nd value of dim_linkCell <= 0"
        if (dim_linkCell(3) <= 0) STOP "ERROR: 3rd value of dim_linkCell <= 0"
    end if
end subroutine check_parsed_calculation

subroutine check_parsed_tasks()
    if ((save_transitory .neqv. .True.) .and. (save_transitory .neqv. .False.)) STOP "ERROR: invalid save_transitory"
    if ((save_observables .neqv. .True.) .and. (save_observables .neqv. .False.)) STOP "ERROR: invalid save_observables"
    if ((save_positions .neqv. .True.) .and. (save_positions .neqv. .False.)) STOP "ERROR: invalid save_positions"
    if ((save_state .neqv. .True.) .and. (save_state .neqv. .False.)) STOP "ERROR: invalid save_state"
    if ((do_pair_correlation .neqv. .True.) .and. (do_pair_correlation .neqv. .False.)) &
        STOP "ERROR: invalid do_pair_correlation"
    if ((do_mean_sqr_displacement .neqv. .True.) .and. (do_mean_sqr_displacement .neqv. .False.)) &
        STOP "ERROR: invalid do_mean_sqr_displacement"
    if ((do_structure_factor .neqv. .True.) .and. (do_structure_factor .neqv. .False.)) &
        STOP "ERROR: invalid do_structure_factor"
end subroutine check_parsed_tasks

subroutine check_parsed_approximation()
    if ((integrator /= 'velocity-Verlet') .and. (integrator /= 'Brownian') .and. (integrator /= 'Monte-Carlo')) then
        STOP "ERROR: unavailable integrator"
    end if

    if (sigma <= 0._pr)     STOP "ERROR: sigma <= 0._pr"
    if (epsilon <= 0._pr)   STOP "ERROR: epsilon <= 0._pr"

    if (integrator == 'Monte-Carlo') then
        if (MC_adjust_step <= 0) STOP "ERROR: MC_adjust_step <= 0"
        if (MC_delta <= 0._pr)   STOP "ERROR: MC_delta <= 0._pr"
    end if

    if  ((integrator == 'velocity-Verlet') .and. (interactions == 'coulomb') .and. (summation /= 'Ewald')) then
        print*, "WARNING: unavailable summation for Coulomb interactions. Switching to 'Ewald' summation"
        summation = 'Ewald'
    end if

    if ((integrator == 'velocity-Verlet') .and. (summation == 'Ewald')) then
        if (kgrid(1) <= 0) STOP "ERROR: 1st value of kgrid <= 0"
        if (kgrid(2) <= 0) STOP "ERROR: 2nd value of kgrid <= 0"
        if (kgrid(3) <= 0) STOP "ERROR: 3rd value of kgrid <= 0"
    end if
end subroutine check_parsed_approximation

subroutine check_parsed_thermostat()
    if ((thermostat_type /= 'rescale') .and. (thermostat_type /= 'Berendsen')) then
        STOP "ERROR: unavailable thermostat"
    end if

    if ((thermostat_type == 'Berendsen').and.(Berendsen_time <= 1e-10)) then
        STOP "ERROR: Berendsen_time <= 1e-10"
    end if
end subroutine check_parsed_thermostat

subroutine check_parsed_MSD()
    if (max_correlation <= 0) STOP "ERROR: max_correlation <= 0"
    if (max_correlation <= real_steps/measuring_jump) print*, "WARNING: max_correlation <= real_steps/measuring_jump"
end subroutine check_parsed_MSD

subroutine check_parsed_pair_correlation()
    ! pair_corr_cutoff must be checked as soon as the periodicity is calculated
    if (pair_corr_bins <= 10) STOP "ERROR: pair_corr_bins <= 10"
end subroutine check_parsed_pair_correlation

END MODULE checkParsingMod