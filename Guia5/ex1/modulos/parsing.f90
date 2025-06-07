MODULE parsing
    use precision
    use subrutinas
    implicit none

    character (len=15)                          :: integrator, type
    character (len=6)                           :: structure
    character (len=12)                          :: summation
    logical                                     :: save_transitory, save_observables
    integer(kind=int_large)                     :: MD_steps, transitory_steps, rescale_steps, dim_linkCell(3)

    ! Namelist blocks
    namelist /physical/ structure, lattice_constant, density, initial_Temp, num_atoms, molar_mass, cell_dim
    namelist /calculation/ MD_steps, transitory_steps, rescale_steps, dt, radius_cutoff, pair_corr_cutoff, pair_corr_bins, summation
    namelist /tasks/ save_transitory, save_observables, do_pair_correlation
    namelist /approximation/ integrator, dim_linkCell, type, sigma, epsilon

    CONTAINS

subroutine set_defaults()
        ! Physical problems' characteristics
        structure           = "random"
        lattice_constant    = 1._pr
        cell_dim            = (/1_int_small,1_int_small,1_int_small/)
        num_atoms           = 0
        initial_Temp        = 1.1_pr
        density             = 0._pr
        molar_mass          = 1._pr

        !Calculation settings
        MD_steps         = 1000
        transitory_steps = 1000
        rescale_steps    = 50
        dt               = 0.005_pr
        radius_cutoff    = 2.5_pr
        pair_corr_cutoff = 4.0_pr
        pair_corr_bins   = 100
        summation        = 'all-vs-all'

        ! Tasks
        save_transitory  = .False.
        save_observables = .False.
        do_pair_correlation = .True.

        ! Potential parameters
        integrator  = 'velocity-Verlet'
        dim_linkCell = (/2,2,2/)
        sigma   = 1._pr
        epsilon = 1._pr

end subroutine set_defaults

subroutine parse_input()
    integer                 :: unit_input
    ! Read from input file
    open(newunit=unit_input, file="input.nml", status="old", action="read")
        read(unit_input, nml=physical)
        read(unit_input, nml=calculation)
        read(unit_input, nml=tasks)
        read(unit_input, nml=approximation)
    close(unit_input)

end subroutine parse_input

END MODULE parsing