!***************************************************************
! Program: ex3.f90
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
! Room for improvement: function used could be optimized for the case x^3 by using f(x) = x*x*x
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex3
    use mtmod,    only: grnd    ! Mersenne -Twister RNG
    use precision
    use formats
    use subrutinas
    use funciones
    use mzranmod
    implicit none

    real(kind=pr)             :: rnd_num, exponent, integral, lower_lim, upper_lim, k
    integer(kind=int_large)   :: i, j, N_iterations, samples, unitnum, max_sample_exponent

    abstract interface
        function funcion(x_x,constants)
            use precision
            implicit none
            real(kind=pr)                :: funcion
            real(kind=pr), intent(in)    :: x_x
            real(kind=pr), intent(in)    :: constants
        end function
    end interface

procedure (funcion), pointer :: function_pointer => null()

!##################################################################################################
!       Input settings
!##################################################################################################
    ! Namelist blocks
    namelist /calculation/ exponent, lower_lim, upper_lim, max_sample_exponent


    ! DEFAULT SETTINGS
        exponent = 3._pr
        lower_lim = 0._pr
        upper_lim = 1._pr
        max_sample_exponent = 1000

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=calculation)
    close(unitnum)

function_pointer => raise

!###################################################################################################
!   Start of the calculations
!###################################################################################################

    open(newunit=unitnum, file="./datos/MC_integral_regular_sampling.out")
    do i = 1, max_sample_exponent
        samples = 10*i
        call MC_grnd_integragtion_regular_sampling(lower_lim, upper_lim, function_pointer, exponent, samples, integral)
        write(unitnum,format_style1) samples, integral, 1._pr/(exponent + 1._pr), abs(1._pr/(exponent + 1._pr) - integral )
    end do
    close(unitnum)

    open(newunit=unitnum, file="./datos/MC_integral_importance_sampling_k2.out")
    k = 2._pr
    do i = 1, max_sample_exponent
        samples = 10*i
        call MC_grnd_integragtion_importance_sampling(lower_lim, upper_lim, function_pointer, exponent, K,samples, integral)
        write(unitnum,format_style1) samples, integral, 1._pr/(exponent + 1._pr), abs(1._pr/(exponent + 1._pr) - integral )
    end do
    close(unitnum)

    open(newunit=unitnum, file="./datos/MC_integral_importance_sampling_k3.out")
    k = 3._pr
    do i = 1, max_sample_exponent
        samples = 10*i
        call MC_grnd_integragtion_importance_sampling(lower_lim, upper_lim, function_pointer, exponent, K,samples, integral)
        write(unitnum,format_style1) samples, integral, 1._pr/(exponent + 1._pr), abs(1._pr/(exponent + 1._pr) - integral )
    end do
    close(unitnum)

end program ex3