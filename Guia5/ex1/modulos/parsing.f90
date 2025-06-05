MODULE parsing
    use precision
    use subrutinas
    implicit none

    character (len=15)                          :: integrator, type
    character (len=6)                           :: structure
    character (len=12)                          :: summation
    logical                                     :: save_transitory, save_observables
    integer(kind=int_huge)                      :: MD_steps, transitory_steps, rescale_steps

    ! Namelist blocks
    namelist /physical/ structure, lattice_constant, density, initial_Temp, num_atoms, molar_mass, cell_dim
    namelist /calculation/ MD_steps, transitory_steps, rescale_steps, dt, radius_cutoff, summation
    namelist /tasks/ save_transitory, save_observables
    namelist /approximation/ integrator, type, sigma, epsilon

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
        summation        = 'all-vs-all'

        ! Tasks
        save_transitory  = .False.
        save_observables = .False.

        ! Potential parameters
        integrator  = 'velocity-Verlet'
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