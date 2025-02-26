MODULE subrutinas
use ISO
use precision
use funciones
implicit none

contains

subroutine three_point_centre(f,x,h,df)
    real (kind=pr)                                          :: f,df
    real (kind=pr), intent (in)                             :: x,h

    df=(f(x+h)-f(x-h))/(2*h)
    !error=h**2*f(eps)/6, optimum h=qbrt(3*eps_machine/M) with (dddf/dddx)<=M

end subroutine three_point_centre

END MODULE