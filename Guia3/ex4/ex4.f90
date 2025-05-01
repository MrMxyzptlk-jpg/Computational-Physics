!***************************************************************
! Program: ex4.f90
! Purpose:
!
!
! Description:
!
! Input: 
!
!
! Output: only if compare_MC is set to true, will the corresponding files be saved in a "./datos" directory with self-explanatory names.
!
!
!
! Room for improvement: function used could be optimized for the case x^3 by using f(x) = x*x*x
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex4
    use precision
    use formats
    use subrutinas
    use funciones
    use mzranmod
    implicit none

    real(pr), allocatable           :: lower_lims(:), upper_lims(:)
    integer(int_huge), allocatable  :: N_points(:)
    integer(int_large)              :: i, dim, unitnum_serial, unitnum_parallel, unitnum
    real (pr)                       :: MC_vol, analytic_vol, trap_vol, MC_stddev
    integer                         :: MC_samples
    logical                         :: do_MC, do_trap, compare_MC
    real(kind=pr)                   :: t_start, t_end, elapsed_time
    character(len=:), allocatable   :: file_MC_serial, file_MC_parallel

    abstract interface
        function funcion(x_x)
            use precision
            implicit none
            real(kind=pr)               ::funcion
            real(kind=pr), intent(in)   :: x_x(:)
        end function
    end interface
    procedure (funcion), pointer :: function_pointer => null()

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ dim, MC_samples, do_MC, do_trap, compare_MC


    ! DEFAULT SETTINGS
    dim = 6     ! Dimension of the N-ball
    MC_samples = 4e7

    do_MC   = .true.
    do_trap = .true.
    compare_MC = .false.

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)

    function_pointer => parametrization
!###################################################################################################
!   Initializing trapezium limits

    if (do_trap) then
        lower_lims  = (/(-1._pr, i = 1, dim-1)/)
        upper_lims  = (/(1._pr, i = 1, dim-1)/)

        select case (dim)
            case (2,3,4,5,7,9,13)
            N_points    = (/(2**(24/(dim-1)),i=1,dim-1)/)
            case default
                do_trap = .false.
        end select
    end if

!###################################################################################################


    if (do_MC) then
        analytic_vol = volume(dim)
        write(*,"(a,i3,a,E16.9)") "Dimensions = ", dim, "  analytic = ", analytic_vol
        print*, "SERIAL RUN"
        call cpu_time(t_start)
        call MC_volume_gaussian_with_uncertainty_serial(dim, MC_samples, function_pointer, MC_vol, MC_stddev)
        call cpu_time(t_end)
        elapsed_time = t_end - t_start
        write(*,"(a,E16.9,a,E16.9,a,E16.9)") "Monte Carlo integration = ", MC_vol, ", Relative error = " &
            , abs(MC_vol-analytic_vol)/analytic_vol, ", Elapsed time = ", elapsed_time

        print*, "PARALLEL RUN"
        t_start = omp_get_wtime()
        call MC_volume_gaussian_with_uncertainty_parallel(dim, MC_samples, function_pointer, MC_vol, MC_stddev)
        t_end = omp_get_wtime()
        elapsed_time = t_end - t_start
        write(*,"(a,E16.9,a,E16.9,a,E16.9,a,E16.9)") "Monte Carlo integration = ", MC_vol, ", Relative error = " &
            , abs(MC_vol-analytic_vol)/analytic_vol, ", Elapsed time = ", elapsed_time, ", Error estimation = ", MC_stddev
    end if

    if (do_trap) then
        analytic_vol = volume(dim)
        write(*,"(a,i3,a,E16.9)") "Dimensions = ", dim, "  analytic = ", analytic_vol
        call cpu_time(t_start)
        trap_vol = integrate_trap(lower_lims, upper_lims, N_points, function_pointer)
        call cpu_time(t_end)
        elapsed_time = t_end - t_start
        write(*,"(a,E16.9,a,E16.9,a,E16.9)") "Traptezium integration  = ", trap_vol, ", Relative error = " &
            , abs(trap_vol-analytic_vol)/analytic_vol, ", Elapsed time = ", elapsed_time
    end if

    if (compare_MC) then
        file_MC_serial = "./datos/MC_serial.out"
        file_MC_parallel = "./datos/MC_parallel.out"

        open(newunit=unitnum_serial, file=file_MC_serial, status="replace")
        open(newunit=unitnum_parallel, file=file_MC_parallel, status="replace")
            write(unitnum_serial,*) "## Volume | Relative error | Calculation time | Estimated error"
            write(unitnum_parallel,*) "## Volume | Relative error | Calculation time | Estimated error"
            do i=2,100
                analytic_vol = volume(i)
                write(*,"(a,i3,a,E16.9)") "Dimensions = ", i, "  analytic = ", analytic_vol

                print*, "SERIAL RUN"
                call cpu_time(t_start)
                call MC_volume_gaussian_with_uncertainty_serial(i, MC_samples, function_pointer, MC_vol, MC_stddev)
                call cpu_time(t_end)
                elapsed_time = t_end - t_start
                write(*,"(a,E16.9,a,E16.9,a,E16.9)") "Monte Carlo integration = ", MC_vol, ", Relative error = " &
                    , abs(MC_vol-analytic_vol)/analytic_vol, ", Elapsed time = ", elapsed_time, ", Error estimation = ", MC_stddev
                write(unitnum_serial,format_style2) MC_vol, abs(MC_vol-analytic_vol)/analytic_vol, elapsed_time &
                , MC_stddev

                print*, "PARALLEL RUN"
                t_start = omp_get_wtime()
                call MC_volume_gaussian_with_uncertainty_parallel(i, MC_samples, function_pointer, MC_vol, MC_stddev)
                t_end = omp_get_wtime()
                elapsed_time = t_end - t_start
                write(*,"(a,E16.9,a,E16.9,a,E16.9)") "Monte Carlo integration = ", MC_vol, ", Relative error = " &
                    , abs(MC_vol-analytic_vol)/analytic_vol, ", Elapsed time = ", elapsed_time, ", Error estimation = ", MC_stddev
                write(unitnum_parallel,format_style2) MC_vol, abs(MC_vol-analytic_vol)/analytic_vol, elapsed_time &
                , MC_stddev

            end do
        close(unitnum_serial)
        close(unitnum_parallel)
    end if

end program ex4
