MODULE subrutinas
use ISO
use precision
use funciones
implicit none

contains

subroutine trapecio(f,a,b,n,integ)
    real (kind=pr), intent (out)                            :: integ
    real (kind=pr), intent (in)                             :: a,b
    real (kind=pr)                                          :: f,h
    integer, intent (in)                                    :: n
    integer                                                 :: ii

    h=(b-a)/real(n,kind=pr)
    integ=0._pr

    do ii=0,n
        if (ii==0.or.ii==n)then
            integ=integ+f(a+real(ii,kind=pr)*h)/2._pr
        else
        integ=integ+f(a+real(ii,kind=pr)*h)
        end if
    end do
    integ=integ*h

    !error(centre)=-h*3*f(eps)12

end subroutine trapecio

subroutine simpson(f,a,b,n,integ)
    real (kind=pr), intent (out)                            :: integ
    real (kind=pr), intent (in)                             :: a,b
    real (kind=pr)                                          :: f,h
    integer, intent (in)                                    :: n
    integer                                                 :: ii

    h=(b-a)/real(n,kind=pr)
    integ=0._pr

    do ii=0,n
        if (ii==0.or.ii==n)then
            integ=integ+f(a+real(ii,kind=pr)*h)
        else if (mod(ii,2)==0.and.ii.ne.0.and.ii.ne.n) then
            integ=integ+2._pr*f(a+real(ii,kind=pr)*h)
        else if (mod(ii,2)==1.and.ii.ne.0.and.ii.ne.n) then
            integ=integ+4._pr*f(a+real(ii,kind=pr)*h)
        end if
    end do
    integ=integ*h/3.0_pr

    !error(centre)=-h**5*f(eps)/90

end subroutine simpson

END MODULE