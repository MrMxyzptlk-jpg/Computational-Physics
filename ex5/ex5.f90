program ex5
    use ISO
    use precision
    use funciones
    use subrutinas
    implicit none

    integer                                     :: i,j,unitnum, unitnum2
    real(kind=pr)                               :: a,b,y0,h, time_euler, time_RK2, time_RK4, count_rate
    real(kind=pr)                               :: euler_avg_prec=0._pr , RK2_avg_prec=0._pr , RK4_avg_prec=0._pr 
    integer(kind=int_huge)                      :: start_euler, start_RK2, start_RK4, finish_euler, finish_RK2, finish_RK4
    real(kind=pr), dimension(:), allocatable    :: exacta,w_euler,t_euler,w_RK2,t_RK2,w_RK4,t_RK4
    character (len=:), dimension(:), allocatable           :: file_out
    character (len=:), allocatable              :: file_times

    abstract interface
    function funcion(x_x,y_y)
        real(kind=8)             ::funcion
        real(kind=8), intent(in) :: x_x,y_y
    end function
    end interface

    procedure (funcion), pointer :: f => null()

    
    f=>f1
    a=0._pr     !lower limit
    b=1._pr     !upper limit
    y0=1._pr    !initial condition
    call SYSTEM_CLOCK(count_rate=count_rate)

    
    file_out = [ "./datos/datos_ex5_1.out", "./datos/datos_ex5_2.out", "./datos/datos_ex5_3.out", "./datos/datos_ex5_4.out", &
    "./datos/datos_ex5_5.out"]
    file_times = "./datos/times.out"
    open(newunit=unitnum2, file=file_times,status="unknown")
    write(unitnum2,"(10x,A15,10x,A10,10x,A13,20x,A8,10x,10x,A13,10x,A8)") "euler precision", "euler time", "RK2 precision", &
    "RK2 time", "RK4 precision", "RK4 time"
    
    do i=1,5
        
        h=10._pr**(-i)    !integration step
        exacta=(/(f2(y0+h*i),i=1,int((b-a)/h))/)    !exact solution

        !Prepare output file
        open(newunit=unitnum, file=file_out(i),status="unknown")
        write(unitnum,"(10x,A1,25x,A1,30x,A5,10x,A14,20x,A3,20x,A14,20x,A3,20x,A14)") "x","y", &
        "Euler", "Error absoluto", "RK2","Error absoluto", "RK4", "Error absoluto"
        !write(unitnum,"(10x,A1,25x,A1,10x,A6,10x,A5,10x,A14,10x,A3,6x,A14,5x,A3,5x,A14)")

        call SYSTEM_CLOCK(start_euler)
        call euler(f,a,b,h,y0,w_euler,t_euler)
        call SYSTEM_CLOCK(finish_euler)
        time_euler = real(finish_euler - start_euler, 8) / real(count_rate, 8)

        call SYSTEM_CLOCK(start_RK2)
        call RK_2(f,a,b,h,y0,w_RK2,t_RK2)
        call SYSTEM_CLOCK(finish_RK2)
        time_RK2 = real(finish_RK2 - start_RK2, 8) / real(count_rate, 8)

        call SYSTEM_CLOCK(start_RK4)
        call RK_4(f,a,b,h,y0,w_RK4,t_RK4)
        call SYSTEM_CLOCK(finish_RK4)
        time_RK4 = real(finish_RK4 - start_RK4, 8) / real(count_rate, 8)

        do j=1,size(t_euler)
            write(unitnum,*) t_euler(j),exacta(j),w_euler(j),abs(w_euler(j)-exacta(j)),&
            w_RK2(j),abs(w_RK2(j)-exacta(j)), w_RK4(j),abs(w_RK4(j)-exacta(j))
            !average precision for each calculation
            euler_avg_prec = euler_avg_prec + abs(w_euler(j)-exacta(j))
            RK2_avg_prec = RK2_avg_prec + abs(w_RK2(j)-exacta(j))
            RK4_avg_prec = RK4_avg_prec + abs(w_RK4(j)-exacta(j))
        end do
        euler_avg_prec = euler_avg_prec /real(size(t_euler),8)
        RK2_avg_prec = RK2_avg_prec /real(size(t_euler),8)
        RK4_avg_prec = RK4_avg_prec /real(size(t_euler),8)

        write(unitnum2,*) euler_avg_prec, time_euler, RK2_avg_prec, time_RK2, RK4_avg_prec, time_RK4

        
        close(unitnum)
        deallocate(w_euler,t_euler,w_RK2,t_RK2,w_RK4,t_RK4, exacta)
    END DO
    deallocate(file_out, file_times)
    close(unitnum2)
    !9 FORMAT("(10x,A1,25x,A1,10x,A6,10x,A14,10x,A4,10x,A3,6x,A14,5x,A4,5x,A3,5x,A14,5x,A4)") !"x","y", &
    !"Euler", "Error absoluto", "time", "RK2","Error absoluto", "time", "RK4", "Error absoluto", "time"

end program ex5