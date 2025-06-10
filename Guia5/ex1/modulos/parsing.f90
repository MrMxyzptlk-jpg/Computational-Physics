MODULE parsing
    use precision
    use subrutinas
    implicit none

    character (len=15)  :: integrator, type
    character (len=6)   :: ensemble
    character (len=12)  :: summation, thermostat_type
    logical             :: save_observables, save_positions, do_structure_factor, do_mean_sqr_displacement
    integer(int_large)  :: MD_steps, transitory_steps, thermostat_steps, dim_linkCell(3), Miller_index(3)
    integer(int_large)  :: correlation_jump

    ! Namelist blocks
    namelist /physical/ structure, lattice_constant, density, initial_Temp_Adim, num_atoms, mass, cell_dim, ensemble
    namelist /calculation/ MD_steps, transitory_steps, thermostat_steps, dt, radius_cutoff, summation, dim_linkCell
    namelist /tasks/ save_transitory, save_positions, save_observables, do_pair_correlation, do_mean_sqr_displacement &
        , do_structure_factor, Miller_index
    namelist /approximation/ integrator, type, sigma, epsilon
    namelist /thermostat/ thermostat_type, Berendsen_time
    namelist /MSD/ max_correlation, correlation_jump
    namelist /pair_correlation/pair_corr_cutoff, pair_corr_bins

    CONTAINS

subroutine set_defaults()
        ! Physical problems' characteristics
        ensemble            = "NVE"
        structure           = "FCC"
        lattice_constant    = 1._pr
        cell_dim            = (/1_int_small,1_int_small,1_int_small/)
        num_atoms           = 0
        initial_Temp_Adim        = 1.1_pr
        density             = 0._pr
        mass                = 1._pr

        !Calculation settings
        MD_steps         = 1000
        transitory_steps = 1000
        thermostat_steps = 50
        dt               = 0.005_pr
        radius_cutoff    = 2.5_pr
        summation        = 'all-vs-all'
        dim_linkCell = (/3,3,3/)

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

        ! Thermostat parameters
        thermostat_type      = 'rescale'
        Berendsen_time  = 0.01_pr

        ! MSD
        max_correlation  = 100
        correlation_jump  = 10

        ! Pair correlation parameters
        pair_corr_cutoff = 4.0_pr
        pair_corr_bins   = 100

end subroutine set_defaults

subroutine parse_input()
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

end subroutine parse_input

END MODULE parsing