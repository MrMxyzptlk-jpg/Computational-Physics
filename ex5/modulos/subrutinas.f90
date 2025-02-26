MODULE subrutinas
use ISO
use precision
use funciones
implicit none

contains

subroutine euler(f,a,b,h,y0,w,t)
    real (kind=pr), dimension(:), allocatable, intent (out) :: w,t
    real (kind=pr), intent (in)                             :: y0,a,b,h
    real (kind=pr)                                          :: f
    integer                                                 :: ii,NN

    NN=int((b-a)/h)

    allocate(w(NN+1))
    t=(/(a+ii*h,ii=0,NN)/)
    w(1)=y0

    do ii=1,NN   
        w(ii+1)=w(ii)+h*f(t(ii),w(ii))   
    end do
    
end subroutine euler

subroutine RK_2(f,a,b,h,y0,w,t)
    real (kind=pr), dimension(:), allocatable, intent (out) :: w,t
    real (kind=pr), intent (in)                             :: y0,a,b,h
    real (kind=pr)                                          :: f
    integer                                                 :: ii,NN

    NN=int((b-a)/h)

    allocate(w(NN+1))
    t=(/(a+ii*h,ii=0,NN)/)
    w(1)=y0

    do ii=1,NN
        w(ii+1)=w(ii)+h*f(t(ii)+h*0.5_pr,w(ii)+h*0.5_pr*f(t(ii),w(ii)))   
    end do
    
end subroutine RK_2

subroutine RK_4(f,a,b,h,y0,w,t)
    real (kind=pr), dimension(:), allocatable, intent (out) :: w,t
    real (kind=pr), intent (in)                             :: y0,a,b,h
    real (kind=pr)                                          :: f,k1,k2,k3,k4
    integer                                                 :: ii,NN

    NN=int((b-a)/h)

    allocate(w(NN+1))
    t=(/(a+ii*h,ii=0,NN)/)
    w(1)=y0

    do ii=1,NN 
        k1=h*f(t(ii),w(ii))
        k2=h*f(t(ii)+h*0.5_pr,w(ii)+0.5_pr*k1)
        k3=h*f(t(ii)+h*0.5_pr,w(ii)+0.5_pr*k2)
        k4=h*f(t(ii+1),w(ii)+k3)
        w(ii+1)=w(ii)+(k1+2._pr*(k2+k3)+k4)/6._pr  
    end do
    
end subroutine RK_4


END MODULE