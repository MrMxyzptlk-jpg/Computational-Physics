program ex3
    use ISO
    use precision
    use funciones
    use subrutinas
    implicit none

    integer                                  :: j ,unitnum
    character(len=:), allocatable            :: file_out
    real(kind=pr)                            :: I, integral_simpson, integral_trapecio
    real(kind=pr), parameter                 :: a=0, b=1
    integer, dimension(:), allocatable       :: n


    !Calculamos el valor exacto de la integral y declaramos los n a utilizar
    I=(f1(b)-f1(a))
    n=(/(2**(j),j=2,8)/)

    !Prepare output file
    file_out="./datos/datos_ex3.out"
    open(newunit=unitnum, file=file_out,status="unknown")
    write(unitnum,"(10x,A1,25x,A1,10x,A6,10x,A24,10x,A14,10x,A23,10x,A14)") "h","n","Exacta","Aproximación (trapecio)", &
    &"Error absoluto", "Aproximación (Simpson)","Error absoluto"

    !nos aseguramos que n sea par, para poder utilizar el método de Simpson 
    do j=1,size(n)
        if (mod(n(j),2)==1) then
            n(j)=n(j)+1
            print*,"n elegido era impar (",n(j)-1,"). Se le sumará 1 para que sea par"
        end if
    end do

    !Loop through values of n, calculate the integral with both approximations and save results to files
    do j=1,size(n)
        call trapecio(f1,a,b,n(j),integral_trapecio)
        call simpson(f1,a,b,n(j),integral_simpson)
        write(unitnum,*)  abs(b-a)/n(j), n(j), I, integral_trapecio,abs(I-integral_trapecio), integral_simpson,&
        &abs(I-integral_simpson)
    end do

    close(unitnum)

    print*, "ex3 done"

    deallocate (file_out)
    !9 FORMAT(15x,A1,22x,A1,25x,A5,15x,A12,10x,A19)
end program ex3