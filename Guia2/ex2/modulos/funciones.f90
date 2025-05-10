MODULE funciones
    USE precision
    USE constantes
    implicit none

    contains

    !funciones de la guía 2

    function analitical_solution(x,t) ! Analitical solution to the non-dimensional problem
        real(kind=pr), dimension(:), intent(in)     :: x
        real(kind=pr), intent(in)                   :: t
        real(kind=pr), dimension(size(x))           :: analitical_solution

        analitical_solution = exp(-t*pi**2) * cos(pi*x)

    end function analitical_solution

    elemental function analitical_solution_test(x,t) ! Analitical solution to the non-dimensional problem
        real(kind=pr), intent(in)               :: t, x
        real(kind=pr)                           :: analitical_solution_test

        analitical_solution_test = exp(-t*pi**2) * cos(pi*x)

    end function analitical_solution_test

    END module