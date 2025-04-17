!***************************************************************
! Program: ex2.f90
! Purpose: numerical solution of the heat equation with null Neumann boundary condition
!
! Description: the input data must be specified in the metric system. The program solves the non-dimensional heat equation and transforms to the proper units using a conversion factor at the end. The jobs are selected with the logical type variables "do_forward_Euler", "do_backward_Euler" and "do_Crank_Nicolson", and all data is written to a directory "./datos" under file names "forward_euler.out", "backward_euler.out" and "crank_nicolson.out" respectively.
!
! Input: 
!
!
! Output: 
!
!
!
! Room for improvement:   
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex2
    use precision
    use funciones
    use subrutinas
    use constantes
    implicit none

    real(kind=pr)                               :: start_time, end_time, dt
    real(kind=pr), dimension(:), allocatable    :: initial_conditions, Y
    real(kind=pr), dimension(:), allocatable    :: conversion_factors
    integer (kind=int_medium)                   :: epochs
    character (len=:), allocatable              :: filename, prefix, file_root_word, suffix
    integer                                     :: i, num_length_points, unitnum
    logical                                     :: do_forward_Euler, do_backward_Euler, do_Crank_Nicolson, do_analitical_sol
    logical                                     :: stability


!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /physical/ conversion_factors
    namelist /calculation/ start_time, end_time, dt, num_length_points, epochs
    namelist /tasks/ do_forward_Euler, do_backward_Euler, do_Crank_Nicolson, do_analitical_sol

    ! DEFAULT SETTINGS
        ! Physical problems' characteristics
        conversion_factors  = (/1._pr, 1._pr/) ! Only relevant information for this model

        !Calculation settings
        start_time          = 0._pr
        end_time            = 1._pr
        dt                  = 0.001_pr
        num_length_points   = 20    ! dx = 1/num_length_points = 0.05
        epochs              = 20   ! Number of iterations before writing to files. Can use int(end_time/dt)  to print only the last iteration but care is needed.

        ! Tasks
        do_forward_Euler    = .False.
        do_backward_Euler   = .False.
        do_Crank_Nicolson   = .False.
        do_analitical_sol   = .False.

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=physical)
        read(unitnum, nml=calculation)
        read(unitnum, nml=tasks)
    close(unitnum)

!##################################################################################################
!      Start of the calculations
!##################################################################################################

    ! Allocate based on input
    allocate(initial_conditions(0:num_length_points))
    initial_conditions = (/(cos(pi * real(i, pr) / real(num_length_points, pr)), i=0, num_length_points)/)

    prefix = "./datos/"
    suffix = ".out"
    start_time = start_time/conversion_factors(2)
    end_time   = end_time/conversion_factors(2)
    dt = dt/conversion_factors(2)

    if (do_forward_Euler) then
        print*, "Calculating Euler Forward"
        call check_stability_criteria(num_length_points,dt, stability)
        if (stability) then
            file_root_word = "forward_euler"
            filename  = prefix//file_root_word//suffix
            call forward_euler_Neumann(initial_conditions, Y, dt, start_time, end_time, filename&
            , conversion_factors, epochs)
            deallocate(filename, file_root_word)
        else
            print*, "Skipping forward Euler calculation"
        end if
    end if

    if (do_backward_Euler) then
        print*, "Calculating Euler Backward"
        file_root_word = "backward_euler"
        filename  = prefix//file_root_word//suffix
        call backward_euler_Neumann(initial_conditions, Y, dt, start_time, end_time, filename&
        , conversion_factors, epochs)
        deallocate(filename, file_root_word)
    end if

    if (do_Crank_Nicolson) then
        print*, "Calculating Crank-Nicolson"
        file_root_word = "crank_nicolson"
        filename  = prefix//file_root_word//suffix
        call crank_nicolson_Neumann(initial_conditions, Y, dt, start_time, end_time, filename&
        , conversion_factors, epochs)
        deallocate(filename, file_root_word)
    end if

    if (do_analitical_sol) then
        print*,"------- Calculating the analitical solution for the final time, t = ", end_time*conversion_factors(2), "seconds"
        file_root_word = prefix//"analitical_sol_end-time_"
        call create_file_name(file_root_word, nint(end_time*conversion_factors(2), int_huge), suffix, filename)

        call analitical_sol(initial_conditions, Y, dt, start_time, end_time, filename, conversion_factors, epochs)
        deallocate(filename, file_root_word)
    end if

end program ex2