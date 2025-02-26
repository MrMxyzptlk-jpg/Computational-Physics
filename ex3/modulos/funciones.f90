MODULE funciones
USE ISO
USE precision
implicit none

contains

!funciones de la guía 1

    function f1(xx)
        real(kind=pr),intent(in)::xx
        real(kind=pr)::f1
        f1=exp(xx)
    end function
END MODULE