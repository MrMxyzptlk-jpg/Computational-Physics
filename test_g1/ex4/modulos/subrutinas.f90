MODULE subrutinas
use, intrinsic :: iso_c_binding
use precision
use constantes
use formats
use funciones
implicit none

contains

subroutine FFTW_forward_write(file_pow,fftw_out,N, factor)
    integer                                                             :: N, i, unit_pow
    real(pr), intent(in)                                                :: factor 
    character(len=:), allocatable, intent(in)                           :: file_pow
    complex(C_double_complex), allocatable, dimension(:), intent(in)    :: fftw_out
    complex(C_double_complex), allocatable, dimension(:)                :: fftw_tmp
    
    fftw_tmp = fftw_out/ real(N, pr)
    open(newunit=unit_pow, file=file_pow, status='replace')
        write(unit_pow, *) "# Frequency | FFTW (Re) | FFTW (Im)"
        do i = -N/2, N/2
            if (i < 0) then
                write(unit_pow, *) factor*real(i, pr), real(conjg(fftw_tmp(-i+1))), &
                aimag(conjg(fftw_out(-i+1)))
            else
                write(unit_pow, *) factor*real(i, pr), real(fftw_tmp(i+1)), aimag(fftw_tmp(i+1))
            end if
        end do
        print*,  '            Values saved in: ', file_pow
        print*, ''
    close(unit_pow)
END SUBROUTINE FFTW_forward_write


subroutine RK_4_array(f,h,y,t, constants)
    real (kind=pr), dimension (:,:), intent (inout)            :: y
    real (kind=pr), dimension (:,:), allocatable               :: k1,k2,k3,k4
    real (kind=pr)                                             :: t
    real (kind=pr), dimension (:), intent (in)                 :: constants
    real (kind=pr), intent (in)                                :: h
    
    interface
    function f(xx,yy, constants)
        use precision
        implicit none
        real(kind=pr), dimension (:,:), allocatable             ::f
        real(kind=pr), intent(in)                               :: xx
        real(kind=pr), dimension (:), intent(in)                :: constants
        real(kind=pr), dimension(:,:), intent(in)               :: yy
    end function
    end interface


    k1=h*f(t,y, constants)
    k2=h*f(t+h*0.5_pr,y+0.5_pr*k1, constants)
    k3=h*f(t+h*0.5_pr,y+0.5_pr*k2, constants)
    k4=h*f(t+h,y+k3, constants)
    y=y+(k1+2._pr*(k2+k3)+k4)/6._pr 
    deallocate(k1,k2,k3,k4)

end subroutine RK_4_array


subroutine DES(f,a,b,h,y0, y, constants)  !Differential Equations System solver 
    real (kind=pr), dimension (:,:), allocatable, intent(out)   :: y
    real (kind=pr), dimension (:,:), intent (in)                :: y0
    real (kind=pr), intent (in)                                 :: a,b,h
    real(kind=pr), dimension (:), intent(in)                    :: constants
    real (kind=pr)                                              :: t
    integer                                                     :: ii,NN
    
    interface
    function f(xx,yy, constants)
        use precision
        implicit none
        real(kind=pr), dimension (:,:), allocatable             ::f
        real(kind=pr), intent(in)                               :: xx
        real(kind=pr), dimension (:), intent(in)                :: constants
        real(kind=pr), dimension(:,:), intent(in)               :: yy
    end function
    end interface

    allocate(y(size(y0,1),size(y0,2)))
    NN=int((b-a)/h)
    t=a
    y=y0

    do ii=1,NN
        call RK_4_array(f,h,y,t, constants)
        t = a + h*ii
    end do

end subroutine DES

subroutine DES_write(file,f,a,b,h,y0, constants)  !Differential Equations System solver 
    real (kind=pr), dimension (:,:), allocatable        :: y
    real (kind=pr), dimension (:,:), intent (in)        :: y0
    real (kind=pr), intent (in)                         :: a,b,h
    real(kind=pr), dimension (:), intent(in)            :: constants
    real (kind=pr)                                      :: t
    integer                                             :: ii,NN,unitnum
    character (len=:), allocatable, intent (in)         :: file
    real(kind=pr), dimension (:), allocatable           :: energy, energy_0, energy_error
    
    interface
    function f(xx,yy, constants)
        use precision
        implicit none
        real(kind=pr), dimension (:,:), allocatable             ::f
        real(kind=pr), intent(in)                               :: xx
        real(kind=pr), dimension (:), intent(in)                :: constants
        real(kind=pr), dimension(:,:), intent(in)               :: yy
    end function
    end interface

    allocate(y(size(y0,1),size(y0,2)))
    NN=int((b-a)/h)
    t=a
    y=y0
    call get_energy(energy_0,y, constants)
    energy = energy_0
    call get_energy_error(energy_0, energy, energy_error)

    ! We write in a file the time, all the first coordinates for all initial conditions, all first derivative of said coordinate
    ! for all initial conditions and so on. This way we utilize arrays of shape MxN to solve a system of N differential equations  
    ! with M different sets of initial conditions.
    open(newunit=unitnum,file=file)
        write(unitnum,*) "## t | y1_1(t) | y1_2(t) | - - - | ym_n-1(t) | ym_n(t) | energy error 1 | - - - | energy error m"
        write(unitnum,format_style) t, y, energy_error  
        do ii=1,NN
            call RK_4_array(f,h,y,t, constants)
            t = a + h*ii
            call get_energy(energy,y, constants)
            call get_energy_error(energy_0, energy, energy_error)
            write(unitnum,format_style) t, y, energy_error  
            !print*, energy_0, energy
        end do
    close(unitnum)

    deallocate(y)


end subroutine DES_write

subroutine Poincare_section(file,f,a,b,h,y0, constants)   
    real (kind=pr), dimension (:,:), allocatable               :: y, y_tmp
    real (kind=pr), dimension (:,:), intent (in)               :: y0
    real (kind=pr), dimension (:), intent (in)                 :: constants
    real (kind=pr), intent (in)                                :: a,b,h
    real (kind=pr)                                             :: t
    integer                                                    :: ii,jj,NN,unitnum
    character (len=:), allocatable, intent (in)                :: file
    
    interface
    function f(xx,yy, constants)
        use precision
        implicit none
        real(kind=pr), dimension (:,:), allocatable             ::f
        real(kind=pr), intent(in)                               :: xx
        real(kind=pr), dimension (:), intent(in)                :: constants
        real(kind=pr), dimension(:,:), intent(in)               :: yy
    end function
    end interface

    allocate(y(size(y0,1),size(y0,2)))
    NN=int((b-a)/h)
    t=a
    y=y0
    y_tmp=y0

    open(newunit=unitnum,file=file, status = 'replace')
        do ii=1,NN
            call RK_4_array(f,h,y,t,constants)
            t = a + h*ii
            do jj = 1, size(y,1)
                if ((y_tmp(jj,3)*y(jj,3) < 0.0_pr) .and. y(jj,4)>0) then
                    write(unitnum,format_style) y(jj,1), y(jj,2)
                end if
            end do
            y_tmp = y
        end do
    close(unitnum)

    deallocate(y)

end subroutine Poincare_section

subroutine get_energy_error(energy_0, energy, energy_error)
    real(kind=pr), dimension(:), intent(in)                 :: energy_0, energy
    real(kind=pr), dimension(:), allocatable, intent(out)   :: energy_error
    logical                                                 :: zero_energy
    integer                                                 :: i

    zero_energy = .False.
    do i = 1, size(energy_0)
        if (energy_0(i) == 0._pr) then
            zero_energy = .True.
            EXIT
        end if
    end do

    if (zero_energy) then 
        energy_error = abs(energy-energy_0)
    else 
        energy_error = abs((energy-energy_0)/energy_0)
    end if

end subroutine get_energy_error

subroutine get_energy(E, y, constants)
    real(kind=pr), dimension (:), intent(in)                :: constants
    real(kind=pr), dimension(:,:), intent(in)               :: y  ! Input: m x 4
    real(kind=pr), dimension (:), allocatable, intent(out)  :: E
    integer                                                 :: i
    
    ! Calculate p2 given q1, p1 and q2 such that E is conserved
    E = sum(y*y,2)/2._pr + constants(1)*y(:,1)*y(:,1)*y(:,3)*y(:,3)
    
    do i = 1, size(E)
        if (E(i) == 0._pr) then
            print*, "Initial energy = 0. Absolute error will be recorded instead of relative error."
            EXIT
        end if
    end do

end subroutine get_energy

subroutine conserve_energy(E, y, constants)
    real(kind=pr), intent(in)                       :: E
    real(kind=pr), dimension (:), intent(in)        :: constants
    real(kind=pr), dimension(:,:), intent(inout)    :: y  ! Input: m x 4

    ! Calculate p2 given q1, p1 and q2 such that E is conserved. Should add a warning if the sqrt argument is negative. This method might not be the best, and a simple regula falsi method could be added after this estimation to ensure the desired conservation of the energy given a certain tolerance.
    y(:,2) = sqrt(2._pr*E - y(:,4)*y(:,4) - y(:,1)*y(:,1) - y(:,3)*y(:,3) - 2_pr*constants(1)*y(:,1)*y(:,1)*y(:,3)*y(:,3))

end subroutine conserve_energy

subroutine generate_initial_conditions(IC, bounds, N)
    integer, intent(in)                                     :: N
    integer                                                 :: ii, jj
    real(kind=pr), dimension(:,:), intent(in)               :: bounds ! 2xm array
    real(kind=pr), dimension(:,:), allocatable, intent(out) :: IC 

    allocate(IC(N, size(bounds,2)))

    do ii = 1, size(bounds,2)
        IC(:,ii) = (/(bounds(1,ii) + (bounds(2,ii) - bounds(1,ii))*real(jj-1,pr)/real(N,pr), jj=1,N)/)
    end do

end subroutine generate_initial_conditions

subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                :: fmt  ! Format descriptor
    character(len=12)                :: x1   ! Temporary string for formatted real
    character(len=:), allocatable    :: prefix, suffix, filename
    real(kind=pr), intent(in)        :: num  ! Input real number

    fmt = '(F8.2)'  ! Format real with 5 decimal places (adjust as needed)
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

END MODULE
