MODULE parsing
    use precision
    use subrutinas
    implicit none

    integer                 :: unit_input
    real(kind=pr)           :: KbT_min, KbT_max, KbT_user, initial_magnetization
    integer(int_large)      :: MC_steps, step_jump, transitory_steps, KbT_steps, lattice_frames
    logical                 :: T_range, save_thermalization, use_absolute_magnetization, do_autocorrelation &
        , save_lattice_evolution, save_transitory, do_binder


    ! Namelist blocks
    namelist /physical/ KbT_min, KbT_max, KbT_steps, T_range, KbT_user, initial_magnetization, x_size, y_size
    namelist /calculation/ MC_steps, step_jump, transitory_steps, use_absolute_magnetization
    namelist /tasks/ save_thermalization, lattice_frames, do_binder , autocorrelation_len_max, do_autocorrelation &
        , save_lattice_evolution, save_transitory

    CONTAINS

subroutine set_defaults()
    ! Physical problems' characteristics
        KbT_user    = 2.2676_pr
        T_range     = .false.
        KbT_min     = 0.1_pr
        KbT_max     = 3.5
        KbT_steps   = 33
        x_size      = 40
        y_size      = 40
        initial_magnetization = 1._pr

    ! Calculation settings
        MC_steps = 1000000
        step_jump = 1
        transitory_steps = MC_steps/2
        use_absolute_magnetization = .true.

    ! Tasks
        save_transitory = .false.
        save_thermalization = .false.
        autocorrelation_len_max = 100
        do_autocorrelation = .false.
        save_lattice_evolution = .false.
        lattice_frames = 100
        do_binder = .false.

end subroutine set_defaults

subroutine parse_input()
    ! Read from input file
    open(newunit=unit_input, file="input.nml", status="old", action="read")
        read(unit_input, nml=physical)
        read(unit_input, nml=calculation)
        read(unit_input, nml=tasks)
    close(unit_input)

end subroutine parse_input

END MODULE parsing