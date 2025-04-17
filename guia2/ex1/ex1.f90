
!***************************************************************
! Program: ex1.f90
! Purpose: numerical solution of the heat equation with null Dirichlet boundary condition
!
! Description: the input data must be specified in the metric system. The program solves the non-dimensional heat equation and transforms to the proper units using a conversion factor at the end. The jobs are selected with the logical type variables "do_forward_Euler", "do_backward_Euler", "do_Crank_Nicolson", "do_analitical_sol" and "do_analitical_conv"
!
! Input: The input data must be specified in an input file named "input.xml" and an example is provided in this repository. All arguments have a default value, and all jobs are set to FALSE as default. 
!
!
! Output: all data is written to a directory "./datos" under file names indicating the task done. 
!
!
!
! Room for improvement: used "initial_condition" as an array with the boundary conditions as well as initial conditions. They should be separated for clarity and is properly implemented in the next program.
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precision
    use funciones
    use subrutinas
    use constantes
    implicit none

    real(kind=pr)                               :: start_time, end_time, dt, initial_Temp, bar_length
    real(kind=pr)                               :: heat_conductivity, specific_heat, density
    real(kind=pr), dimension(:), allocatable    :: initial_conditions, Y
    real(kind=pr), dimension(:), allocatable    :: conversion_factors
    real(kind=pr)                               :: t_start, t_end, elapsed_time
    character (len=:), allocatable              :: filename, prefix, file_root_word, suffix
    integer(kind=int_huge), dimension(:), allocatable   :: max_analitical_terms
    integer(kind=int_huge)                      :: num_analitical_terms
    integer(kind=int_huge)                      :: epochs
    integer                                     :: i, num_length_points, unitnum
    logical                                     :: do_forward_Euler, do_backward_Euler, do_Crank_Nicolson, do_analitical_sol
    logical                                     :: do_analitical_conv, do_conv_initial_conditions, do_conv_end_time
    logical                                     :: stability

!##################################################################################################
!       Input settings
!##################################################################################################

    ! Namelist blocks
    namelist /physical/ bar_length, heat_conductivity, specific_heat, density, initial_Temp
    namelist /calculation/ start_time, end_time, dt, num_length_points, epochs
    namelist /tasks/ do_forward_Euler, do_backward_Euler, do_Crank_Nicolson, do_analitical_sol, num_analitical_terms &
        , do_analitical_conv, do_conv_initial_conditions,  do_conv_end_time

    ! DEFAULT SETTINGS
        ! Physical problems' characteristics
        bar_length          = 1._pr     ! Will be used to re-dimensionalize the problem
        initial_Temp        = 100._pr
        heat_conductivity   = 237._pr
        specific_heat       = 900._pr
        density             = 2700._pr

        !Calculation settings
        start_time          = 0._pr
        end_time            = 1800._pr
        dt                  = 0.3_pr
        num_length_points   = 100       ! dx = 1/(num_length_points + 1)
        epochs              = 300       ! Number of iterations before writing to files. Can use int(end_time/dt)  to print only the last iteration but care is needed.

        ! Tasks
        do_forward_Euler    = .False.
        do_backward_Euler   = .False.
        do_Crank_Nicolson   = .False.
        do_analitical_sol   = .False.
            num_analitical_terms= 1e2
        do_analitical_conv  = .False.
            do_conv_initial_conditions = .False.
            do_conv_end_time = .False.
            max_analitical_terms= (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100/)  ! NOT READ FROM FILE


    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=physical)
        read(unitnum, nml=calculation)
        read(unitnum, nml=tasks)
    close(unitnum)

!##################################################################################################
!      Start of the calculations
!##################################################################################################

    ! Set files' name defaults, initial conditions and adimensionalize variables
    prefix = "./datos/"
    suffix = ".out"
    initial_conditions = (/0._pr,(initial_Temp, i=1, num_length_points-1),0._pr/)
    conversion_factors = (/bar_length, (specific_heat*density*bar_length*bar_length)/heat_conductivity /)
    start_time = start_time/conversion_factors(2)
    end_time   = end_time/conversion_factors(2)
    dt = dt/conversion_factors(2)

    if (do_forward_Euler) then
        print*,"-------------------- Calculating with forward Euler --------------------"
        call check_stability_criteria(num_length_points,dt, stability)
        if (stability) then
            call cpu_time(t_start)

            file_root_word = "forward_euler"
            filename  = prefix//file_root_word//suffix
            call forward_euler(initial_conditions, Y, dt, start_time, end_time, filename, conversion_factors, epochs)
            deallocate(filename, file_root_word)

            call cpu_time(t_end)
            elapsed_time = t_end - t_start
            print*, "Time elapsed:", elapsed_time, "seconds"
        else
            print*, "Skipping forward Euler calculation"
        end if
    end if

    if (do_backward_Euler) then
        call cpu_time(t_start)

        print*,"-------------------- Calculating with backward Euler --------------------"
        file_root_word = "backward_euler"
        filename  = prefix//file_root_word//suffix
        call backward_euler(initial_conditions, Y, dt, start_time, end_time, filename, conversion_factors,epochs)
        deallocate(filename, file_root_word)

        call cpu_time(t_end)
        elapsed_time = t_end - t_start
        print*, "Time elapsed:", elapsed_time, "seconds"
    end if

    if (do_Crank_Nicolson) then
        call cpu_time(t_start)

        print*,"-------------------- Calculating with Crank-Nicolson --------------------"
        file_root_word = "crank_nicolson"
        filename  = prefix//file_root_word//suffix
        call crank_nicolson(initial_conditions, Y, dt, start_time, end_time, filename, conversion_factors,epochs)
        deallocate(filename, file_root_word)

        call cpu_time(t_end)
        elapsed_time = t_end - t_start
        print*, "Time elapsed:", elapsed_time, "seconds"
    end if

    if (do_analitical_sol) then
        print*,"------- Calculating the analitical solution to the sum term N=", num_analitical_terms, " and time = "&
        , end_time*conversion_factors(2), "seconds"
        file_root_word = prefix//"analitical_sol_end-time_"
        call create_file_name(file_root_word, nint(end_time*conversion_factors(2), int_huge), suffix, filename)

        call analitical_sol(initial_conditions, Y, end_time, filename, conversion_factors,num_analitical_terms, initial_Temp)
        deallocate(filename, file_root_word)
    end if

    if (do_analitical_conv) then
        if (do_conv_initial_conditions) then
            do i = 1 , size(max_analitical_terms)
                print*,"------- Calculating the analitical solution to the sum term N=", max_analitical_terms(i), "--------"
                file_root_word = prefix//"convergence/analitical_conv_start-time_0_N_"
                call create_file_name(file_root_word, max_analitical_terms(i), suffix, filename)

                call converge_analitical_sol(initial_conditions, Y, start_time, filename, conversion_factors&
                ,max_analitical_terms(i), initial_Temp)
                deallocate(filename, file_root_word)
            end do
        end if
        if (do_conv_end_time) then
            do i = 1 , size(max_analitical_terms)
                ! Create appropriate file name with end time and maximum sum term used
                print*,"Calculating the analitical solution to the sum term N=", max_analitical_terms(i), " and time = "&
                , end_time*conversion_factors(2), "seconds"
                file_root_word = prefix//"convergence/analitical_conv_end-time_"
                deallocate(suffix)
                suffix = "_N_"
                call create_file_name(file_root_word, nint(end_time*conversion_factors(2), int_huge), suffix, filename)
                deallocate(file_root_word, suffix)
                suffix = ".out"
                file_root_word = filename
                deallocate(filename)
                call create_file_name(file_root_word, max_analitical_terms(i), suffix, filename)

                call analitical_sol(initial_conditions, Y, end_time, filename, conversion_factors, max_analitical_terms(i)&
                , initial_Temp)
                deallocate(filename, file_root_word)
            end do
        end if
    end if


end program ex1
