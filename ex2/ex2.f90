program ex2
    use ISO
    use precision
    use funciones
    use subrutinas
    implicit none
    real(kind=pr)                            :: x,h,df
    real(kind=pr), dimension(:), allocatable :: a
    integer                                  :: i,unitnum
    character(len=:), allocatable            ::file_ex1

    !guardamos los h que utilizaremos en un vector y declaramos el punto en el que evaluaremos la derivada
    a=(/(10._pr**(-i),i=1,20)/)
    x=1._pr

    file_ex1="./datos/datos_ex1.out"
    open(newunit=unitnum, file=file_ex1,status="unknown")
    write(unitnum,"(3x,A1,10x,A1,19x,A5,9x,A11,10x,A19)") "h","x","df(x)","aproximación","Error absoluto en y"
    print*,"           h                        x                        df(x)                aproximación              Error abs&
    &oluto en y"
    do i=1,size(a)
        h=a(i)
        call three_point_centre(f1,x,h,df)
        write(unitnum,*) h,x,df1(x),df,abs(df-df1(x))
        print*,h,x,df1(x),df,abs(df-df1(x))
    end do
    close(unitnum)

    h=(3._pr*1.E-15_pr/ddf1(2.1_pr))**(1._pr/3._pr)
    call three_point_centre(f1,x,h,df)
    Print*, "h optimo:"
    print*,h,x,df1(x),df,abs(df-df1(x))

    deallocate (file_ex1)
    call system ("gnuplot -p plot.gp")

end program ex2