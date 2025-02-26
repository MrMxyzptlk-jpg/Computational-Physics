program ex4
    use ISO
    use precision
    use constantes, only: pi
    use funciones
    use subrutinas
    implicit none

    real (kind=pr), dimension (:), allocatable   :: x, w
    real (kind=pr)                               :: a = 0._pr,b = 1._pr
    real (kind=pr)                               :: integral, integral_gauss,integral_trapecio, integral_simpson
    character (len=:), allocatable               :: file_out
    integer(kind=int_large)                      :: j, unitnum
    integer, dimension(:), allocatable           :: N_points


    abstract interface
    function funcion(x_x,y_y)
        real(kind=8), dimension (:), allocatable      ::funcion
        real(kind=8), intent(in)                 :: x_x
        real(kind=8), dimension(:), intent(in)   :: y_y
    end function  

    end interface

    procedure (funcion), pointer :: f => null()
    procedure (funcion), pointer :: g => null()
    procedure (funcion), pointer :: z => null()
    

    !Calculamos el valor exacto de la integral y declaramos los n a utilizar
    integral=-(f1(b)-f1(a))
    N_points=(/2,(10*2**(j),j=0,8)/)


    !Prepare output file
    file_out="./datos/datos_ex4.out"
    open(newunit=unitnum, file=file_out,status="unknown")
    write(unitnum,"(10x,A1,25x,A1,10x,A6,10x,A24,10x,A14,10x,A23,6x,A14,5x,A30,5x,A14)") "h","n","Exacta", &
    "Aproximación (trapecio)", "Error relativo", "Aproximación (Simpson)","Error relativo", "Aproximación (Gauss Legendre)", &
    "Error relativo"

    !nos aseguramos que n sea par, para poder utilizar el método de Simpson y el resto 
    do j=1,size(N_points)
        if (mod(N_points(j),2)==1) then
            N_points(j)=N_points(j)+1
            print*,"n elegido era impar (",N_points(j)-1,"). Se le sumará 1 para que sea par"
        end if
    end do

    !Loop through values of n, calculate the integral with both approximations and save results to files
    do j=1,size(N_points)

        call trapecio(f1,a,b,N_points(j),integral_trapecio)

        call simpson(f1,a,b,N_points(j),integral_simpson)

        integral_gauss = 0._pr
        allocate(x(N_points(j)),w(N_points(j)))
        call gauss_integral(N_points(j), 0, a, b, x, w, f1, integral_gauss)

        write(unitnum,*)  abs(b-a)/N_points(j), N_points(j), integral, integral_trapecio, abs(integral-integral_trapecio)/integral,&
        &integral_simpson, abs(integral-integral_simpson)/integral, integral_gauss, abs(integral-integral_gauss)/integral
        deallocate(x,w)
    end do

    close(unitnum)

    print*, "ex4 done"

    deallocate (file_out)
    !9 FORMAT(15x,A1,22x,A1,25x,A5,15x,A12,10x,A19)


end program ex4