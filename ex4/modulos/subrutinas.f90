MODULE subrutinas
use ISO
use precision
use funciones
implicit none

contains

subroutine gauss_quadrature(npts, job, a, b, x, w)
    implicit none
    integer, intent(in)     :: npts,job
    real (pr), intent(in)   :: a,b
    real (pr), intent(out)  :: x(0:npts),w(0:npts)
    real (pr)               :: t,t1,pp,p1,p2,p3,eps
    real (pr)               :: pi
    real (pr)               :: xi
    integer                 :: m,i,j
  !
    pi = acos(-1._pr)
    m = 0; i = 0; j = 0
    t = 0._pr;  t1 = 0._pr; pp = 0._pr; 
    p1 = 0._pr; p2 = 0._pr; p3 = 0._pr
    eps = 3.e-14_pr                ! Accuracy: ********ADJUST THIS*********!
    m = (npts + 1)/2
  !       do i=1,m + 1
    do i=1,m + 1
  !          t = cos(pi*(dble(i) - 0.25_pr)/(dble(npts) + 0.5_pr) )
       t = cos(pi*(dble(i) - 0.25_pr)/(dble(npts+1) + 0.5_pr) )
       t1 = 1._pr
       do while( (abs(t - t1) ) >= eps)
          p1 = 1._pr ; p2 = 0._pr
  !
          do j=1,npts + 1
             p3 = p2; p2 = p1
  !                p1 = ( (2._pr*dble(j) - 1._pr)*t*p2 - (dble(j) - 1._pr)*p3)/(dble(j) )
             p1 = ( (2._pr*dble(j) - 1._pr)*t*p2 - dble(j-1)*p3)/(dble(j) )
          enddo
  !
  !             pp = dble(npts)*(t*p1 - p2)/(t*t - 1._pr)
          pp = dble(npts+1)*(t*p1 - p2)/(t*t - 1._pr)
          t1 = t
          t = t1 - p1/pp
       enddo
       x(i - 1) = - t; 
  !          x(npts - i) = t
       x(npts + 1 - i) = t
       w(i - 1) = 2._pr/( (1._pr - t*t)*pp*pp)
  !          w(npts - i) = w(i - 1)
       w(npts + 1 - i) = w(i - 1)
  !              print *," x(i - 1)", x(i - 1) , " w " , w(npts - i),i
    enddo
    if (job == 0) then  !Standard Gauss-Legendre Quadrature
       do i=0, npts
          x(i) = x(i)*(b - a)/2._pr + (b + a)/2._pr
          w(i) = w(i)*(b - a)/2._pr
       enddo
    endif
    if (job == 1) then  !Gauss-Radau Quadrature
       do i=0, npts
          xi = x(i)
          x(i) = a*b*(1._pr + xi) / (b + a - (b - a)*xi)
          w(i) = w(i)*2._pr*a*b*b/( (b + a - (b - a)*xi)*(b + a - (b - a)*xi) )
       enddo
    endif
    if (job == 2) then  !Gauss-Lobatto Quadrature
       do i=0, npts
          xi = x(i)
          x(i) = (b*xi + b + a + a) / (1._pr - xi)
          w(i) = w(i)*2.*(a + b)/( (1._pr - xi)*(1. - xi) )
       enddo
    endif
  !
    return
End subroutine gauss_quadrature

subroutine gauss_integral(npts, job, a, b, x, w, f, integral_gauss)
    implicit none
    integer, intent(in)     :: npts,job
    real (kind=pr)          :: f
    real (pr), intent(in)   :: a,b
    real (pr), intent(out)  :: x(0:npts),w(0:npts), integral_gauss
    integer                 :: i

    ! Compute Gauss-Legendre nodes and weights
    call gauss_quadrature(npts, job, a, b, x, w)

    ! Perform the weighted sum
    do i = 0, npts
        integral_gauss = integral_gauss + w(i) * f(x(i))
    end do

End subroutine gauss_integral

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