MODULE parsing
    use precision
    use subrutinas
    implicit none

    character (len=15)                          :: type
    character (len=6)                           :: structure
    logical                                     :: do_velocity_verlet
    real(kind=pr)                               :: sigma, epsilon
    real(kind=pr)                               :: start_time, end_time, dt
    real(kind=pr)                               :: radius_cutoff, initial_Temp, density, lattice_constant, molar_mass
    integer(kind=int_huge)                      :: epochs

    ! Namelist blocks
    namelist /physical/ structure, lattice_constant, density, initial_Temp, num_atoms, molar_mass, cell_dim
    namelist /calculation/ start_time, end_time, dt, radius_cutoff, epochs
    namelist /tasks/ do_velocity_verlet
    namelist /approximation/ type, sigma, epsilon

    CONTAINS

subroutine set_defaults()
        ! Physical problems' characteristics
        structure           = "random"
        lattice_constant    = 1._pr
        cell_dim            = (/1_int_small,1_int_small,1_int_small/)
        num_atoms           = 0
        initial_Temp        = 1.1_pr
        density             = 0.8_pr
        molar_mass          = 1._pr

        !Calculation settings
        start_time      = 0._pr
        end_time        = 1800._pr
        dt              = 0.3_pr
        radius_cutoff   = 2.5_pr
        epochs          = 300

        ! Tasks
        do_velocity_verlet  = .True.

        ! Potential parameters
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