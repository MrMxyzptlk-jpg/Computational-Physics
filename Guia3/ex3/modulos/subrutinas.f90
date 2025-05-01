MODULE subrutinas
    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int32, real64
    USE precision
    use mzranmod
    implicit none

    contains

    subroutine MC_grnd_integragtion_regular_sampling(lower_lim, upper_lim, f, exponent, N_samples, integral)
        real (kind=pr), intent (in)                   :: lower_lim, upper_lim, exponent
        real (kind=pr), intent (out)                  :: integral
        integer (kind = int_large), intent (in)       :: N_samples
        integer (kind = int_large)                    :: i
        real (kind=pr)                                :: x, rnd_num

        interface
            function f(xx, constants)
                use precision
                real(kind=pr)                :: f
                real(kind=pr), intent(in)    :: xx
                real(kind=pr), intent(in)                :: constants
            end function
        end interface

        integral = 0._pr
        do i = 1, N_samples
            rnd_num = rmzran()
            x = lower_lim + rnd_num*(upper_lim - lower_lim)
            integral = integral + f(x, exponent)
        end do
        integral = integral*(upper_lim - lower_lim)/real(N_samples,pr)

    end subroutine MC_grnd_integragtion_regular_sampling

    subroutine MC_grnd_integragtion_importance_sampling(lower_lim, upper_lim, f, exponent, K, N_samples, integral)
        real (kind=pr), intent (in)                   :: lower_lim, upper_lim, exponent, K
        real (kind=pr), intent (out)                  :: integral
        integer (kind = int_large), intent (in)       :: N_samples
        integer (kind = int_large)                    :: i
        real (kind=pr)                                :: x, rnd_num

        interface
            function f(xx, constants)
                use precision
                real(kind=pr)                :: f
                real(kind=pr), intent(in)    :: xx
                real(kind=pr), intent(in)                :: constants
            end function
        end interface

        integral = 0._pr
        do i = 1, N_samples
            rnd_num = rmzran()
            x = lower_lim + rnd_num*(upper_lim - lower_lim)
            x = x**(1._pr/(K+1._pr))
            integral = integral + f(x, exponent)/((k + 1._pr)*x**k)
        end do
        integral = integral*(upper_lim - lower_lim)/real(N_samples,pr)

    end subroutine MC_grnd_integragtion_importance_sampling

END module subrutinas
