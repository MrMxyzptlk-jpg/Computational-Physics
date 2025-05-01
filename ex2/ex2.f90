!***************************************************************
! Program: ex2.f90
! Purpose:
!
!
! Description:
!
! Input:
!
!
! Output:
!
!
!
! Room for improvement: the subroutines can be vectoried and parallelized for farter calculations, given that every random walk is independent from the rest.
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex2
    use mtmod,    only: grnd    ! Mersenne -Twister RNG
    use precision
    use formats
    use mzranmod
    use funciones
    use subrutinas
    implicit none

    real(kind=pr)                                       :: rnd_num, mean_quadratic_displacement
    integer(kind=int_large), dimension(:), allocatable  :: lattice_position
    real(kind=pr), dimension(:), allocatable            :: cuadrants
    integer(kind=int_small)                             :: dim, step, k, min_steps, max_steps
    integer(kind=int_large)                             :: i, j, N_iterations, n, unitnum, user_seed
    logical                                             :: do_borders, do_msd
    character(len=4)                                    :: method
    character(len=:), allocatable                       :: file_name


!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ dim, N_iterations, method, do_borders, do_msd,  user_seed, min_steps, max_steps

    ! DEFAULT SETTINGS
        N_iterations = 100000
        dim          = 2
        method       = "MZ" ! method = MT (Mersenne - Twister), MZ ( Marsaglia - Zaman), ran2 (from Numerical Recipes)
        do_msd       = .true.
        do_borders   = .true.
        user_seed    = -123
        min_steps    = 1
        max_steps    = 100

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)

!##################################################################################################
!      Start of the calculations
!##################################################################################################

    ! Set up RNG with selected (method and seed if specified)
    call set_rng(method, user_seed)


    if (do_msd) then
        print*, "----------- Calculating Mean Cuadratic Displacement with ", method, " ---------"
        print*, ""
        call create_file_name("./datos/mean_quadratic_displacement_",dim,"d_"//trim(method)//".out", file_name)
        open(newunit=unitnum, file=file_name)
        do step = 1, 100
            call get_mean_square_displacement(N_iterations, dim, step, mean_quadratic_displacement)
            write(unitnum,format_style1) step, mean_quadratic_displacement
        end do
        close(unitnum)
        deallocate(file_name)
    end if

!###########################################################################################################
!   Including cuadrants' borders:
!###########################################################################################################


    if (do_borders) then
        print*, "----------- Calculating Cuadrants Probabilities with ", method, " ---------"
        print*, ""
        call create_file_name("./datos/cuadrants_",dim,"d_"//trim(method)//".out", file_name)
        open(newunit=unitnum, file=file_name)
        do step = 1, 100
            call get_cuadrant_freq(N_iterations, dim, step, cuadrants)
            write(unitnum,format_style1) step, cuadrants
        end do
        close(unitnum)
        deallocate(file_name)
    end if

end program ex2