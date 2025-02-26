MODULE funciones
USE ISO
USE prec_num
USE constantes, only : pi
implicit none

contains

!funciones de la guía 1

    function f1(xx)
        real(kind=wp),intent(in)::xx
        real(kind=wp)::f1
        f1=exp(-xx)
    end function

    function f2(xx)
        real(kind=wp),intent(in)::xx
        real(kind=wp)::f2
        f2=exp(-xx)
    end function

    function f3(xx,yy)
        real(kind=wp), intent(in)                        ::xx
        real(kind=wp),dimension(:), intent(in)           ::yy
        real(kind=wp), dimension(:), allocatable         ::f3

        allocate(f3(size(yy)))
        f3(1)=-10._wp*sin(yy(2))+xx*0
        f3(2)=yy(1)
    end function

    function f4(xx,yy)
        real(kind=wp),intent(in)                          ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::f4

        allocate(f4(size(yy)))
        f4(1)=-10._wp*yy(2)+xx*0
        f4(2)=yy(1)
    end function

    function f5(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::f5

        allocate(f5(size(yy)))
        f5(1)=-0.5_wp*yy(1)*yy(2)/1000._wp+0*xx
        f5(2)=0.5_wp*yy(1)*yy(2)/1000._wp-0.1_wp*yy(2)
        f5(3)=0.1_wp*yy(2)
    end function

    function f6(xx,yy)
        real(kind=wp),intent(in)::xx,yy
        real(kind=wp)::f6
        f6=2._wp*yy*(1._wp-yy/200._wp)+0*xx
    end function

    function f7(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::f7

        allocate(f7(size(yy)))
        f7(1)=-yy(1)**2
        f7(2)=yy(1)
    end function

    function f8(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::f8

        allocate(f8(size(yy)))
        f8(1)=2._wp*yy(2)**3
        f8(2)=yy(1)
    end function

    function f9(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::f9

        allocate(f9(size(yy)))
        f9(1)=((xx*yy(1))**2-9._wp*yy(2)**2+4._wp*xx**6)/(xx**5)
        f9(2)=yy(1)
    end function

!función auxiliar 

    function z7(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp), dimension(:), intent(in)           ::yy
        real(kind=wp), dimension(:), allocatable          ::z7

        allocate(z7(size(yy)))
        z7(1)=-2._wp*yy(3)*yy(1)
        z7(2)=yy(1)
        z7(3)=-yy(3)**2
        z7(4)=yy(3)
    end function

    function z8(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::z8

        allocate(z8(size(yy)))
        z8(1)=6._wp*yy(3)**2*yy(2)
        z8(2)=yy(1)
        z8(3)=2._wp*yy(4)**3
        z8(4)=yy(3)
    end function

    function z9(xx,yy)
        real(kind=wp), intent(in)                         ::xx
        real(kind=wp),dimension(:), intent(in)            ::yy
        real(kind=wp), dimension(:), allocatable          ::z9

        allocate(z9(size(yy)))
        z9(1)=-9._wp*yy(4)*yy(2)/(xx**5)+yy(3)*2._wp/(xx**3)*yy(1)
        z9(2)=yy(1)
        z9(3)=((xx*yy(3))**2-9._wp*yy(4)**2+4._wp*xx**6)/(xx**5)
        z9(4)=yy(3)

    end function

!función real 

    function g1(xx)
        real(kind=wp),intent(in)::xx
        real(kind=wp)::g1
        g1=(1._wp+2._wp*pi/(1._wp+4._wp*pi**2))*exp(-xx)+(sin(2._wp*pi*xx)-2._wp*pi*cos(2._wp*pi*xx))/(1._wp+4._wp*pi**2)
    end function

    function g4 (xx)
        real(kind=wp),intent(in)::xx
        real(kind=wp)::g4
        g4=cos(sqrt(10._wp)*xx)
    end function

END module