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
    

subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                :: fmt  ! Format descriptor
    character(len=12)                :: x1   ! Temporary string for formatted real
    character(len=:), allocatable    :: prefix, suffix, filename
    real(kind=pr), intent(in)        :: num  ! Input real number

    fmt = '(F8.5)'  ! Format real with 5 decimal places (adjust as needed)
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

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

subroutine DES_write(file,f,a,b,h,y0, constants)  !Differential Equations System solver 
    real (kind=pr), dimension (:,:), allocatable               :: y
    real (kind=pr), dimension (:,:), intent (in)               :: y0
    real (kind=pr), intent (in)                                :: a,b,h
    real(kind=pr), dimension (:), intent(in)                    :: constants
    real (kind=pr)                                             :: t
    integer                                                    :: ii,NN,unitnum
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

    ! We write in a file the time, all the first coordinates for all initial conditions, all first derivative of said coordinate
    ! for all initial conditions and so on. This way we utilize arrays of shape MxN to solve a system of N differential equations  
    ! with M different sets of initial conditions.
    open(newunit=unitnum,file=file)
        write(unitnum,format_style) t,y     ! writing "t | y1_1(t) | y1_2(t) | - - - | ym_n-1(t) | ym_n(t) |".
        do ii=1,NN
            call RK_4_array(f,h,y,t, constants)
            t = a + h*ii
            write(unitnum,format_style) t,y 
        end do
    close(unitnum)

    deallocate(y)


end subroutine DES_write


subroutine get_cart_coordinates(file,times,angular_trajectories, cartesian_trajectories, constants, initial_condition_index)
    real (kind=pr), dimension (:,:), allocatable, intent(out)  :: angular_trajectories, cartesian_trajectories
    real (kind=pr), dimension (:), allocatable, intent(out)    :: times
    real (kind=pr), dimension (:), intent (in)                 :: constants
    integer (kind=int_large)                                   :: i
    integer                                                    :: initial_condition_index
    integer                                                    :: io, lines, unitnum
    character (len=:), allocatable, intent (in)                :: file
    
    open (newunit=unitnum, file=file, status='old', access='sequential', form='formatted', action='read' )
        lines = 0 
        DO
            read(unitnum,*,iostat=io)
            if (io/=0) EXIT
            lines = lines + 1
        END DO
        allocate(angular_trajectories(lines,8), cartesian_trajectories(lines,4), times(lines))
        
        rewind(unitnum)
        DO i =1,lines
            read(unitnum,format_style) times(i),angular_trajectories(i,:)
        END DO
        cartesian_trajectories(:,1) = sin(angular_trajectories(:,initial_condition_index))   ! x1
        cartesian_trajectories(:,2) = -cos(angular_trajectories(:,initial_condition_index)) ! y1
        cartesian_trajectories(:,3) = sin(angular_trajectories(:,initial_condition_index)) + &
            constants(2)*sin(angular_trajectories(:,initial_condition_index + 4))  ! x2
        cartesian_trajectories(:,4) = -cos(angular_trajectories(:,initial_condition_index)) - &
            constants(2)*cos(angular_trajectories(:,initial_condition_index + 4))  ! y2
    close(unitnum)

END subroutine get_cart_coordinates

subroutine get_ang_coordinates(file,times,angular_trajectories)
    real (kind=pr), dimension (:,:), allocatable, intent(out)  :: angular_trajectories
    real (kind=pr), dimension (:), allocatable, intent(out)    :: times
    integer (kind=int_large)                                   :: i
    integer                                                    :: io, lines, unitnum
    character (len=:), allocatable, intent (in)                :: file
    
    open (newunit=unitnum, file=file, status='old', access='sequential', form='formatted', action='read' )
        lines = 0 
        DO
            read(unitnum,*,iostat=io)
            if (io/=0) EXIT
            lines = lines + 1
        END DO
        allocate(angular_trajectories(lines,8), times(lines))
        
        rewind(unitnum)
        DO i =1,lines
            read(unitnum,format_style) times(i),angular_trajectories(i,:)
        END DO
    close(unitnum)

END subroutine get_ang_coordinates

subroutine get_generalized_momenta(y,constants,momenta)  ! Used to make the code clearer
    real (kind=pr), dimension (:,:), allocatable, intent (in)   :: y
    real (kind=pr), dimension (:,:)                             :: momenta
    real (kind=pr), dimension (:), intent (in)                  :: constants
    
    momenta(:,1) = compute_p1(y,constants)
    momenta(:,2) = compute_p2(y,constants)

END subroutine get_generalized_momenta

subroutine Poincare_section(file,f,a,b,h,y0, constants)   
    real (kind=pr), dimension (:,:), allocatable               :: y, y_tmp, momenta
    real (kind=pr), dimension (:), allocatable                 :: energies
    real (kind=pr), dimension (:,:), intent (in)               :: y0
    real (kind=pr), dimension (:), intent (in)                 :: constants
    real (kind=pr), intent (in)                                :: a,b,h
    real (kind=pr)                                             :: t, angle
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
    allocate(momenta(size(y,1),2))
    NN=int((b-a)/h)
    t=a
    y=y0
    y_tmp=y0

    open(newunit=unitnum,file=file, status = 'replace')
        do ii=1,NN
            call RK_4_array(f,h,y,t,constants)
            t = a + h*ii
            y(:,1) = mod(y(:,1)+pi, 2._pr*pi) - pi
            y(:,1) = mod(y(:,1)-pi, 2._pr*pi) + pi
            y(:,3) = mod(y(:,3)+pi, 2._pr*pi) - pi
            y(:,3) = mod(y(:,3)-pi, 2._pr*pi) + pi
            call get_generalized_momenta(y,constants,momenta)
            do jj = 1, size(y,1)
                if ((y_tmp(jj,1)*y(jj,1) <= 0.0_pr) .and. momenta(jj,1)>0) then
                !    angle = mod(y(jj,3)+pi, 2._pr*pi) - pi
                !    angle = mod(angle-pi, 2._pr*pi) + pi
                    write(unitnum,format_style) y(jj,3), momenta(jj,2)
                end if
            end do
            y_tmp = y
        end do
    close(unitnum)

    deallocate(y, momenta)

end subroutine Poincare_section

subroutine conserve_energy(E, y, constants)
    real(kind=pr), intent(in)                        :: E
    real(kind=pr), dimension (:), intent(in)          :: constants
    real (kind=pr), dimension (:), allocatable       :: p1, p2
    real(kind=pr), dimension(:,:), intent(inout)     :: y  ! Input: m x 4
    
    p2 = compute_p2(y,constants) 
    ! Calculate p1 given q1, p2 and q2 such that E is conserved
    p1 = p2*cos(y(:,1)-y(:,3))/constants(2) + sqrt((2._pr + constants(1) - constants(1)*cos(2._pr*(y(:,1) - y(:,3)))) &
    *(2._pr*constants(2)*constants(2)*constants(1)*(E + constants(3)*((1._pr + constants(1))*cos(y(:,1))+constants(2)*constants(1) &
    *cos(y(:,3)))) - p2*p2)) /(constants(2)*sqrt(constants(1)*2._pr))

    y(:,2) = (constants(2)*p1 - p2*cos(y(:,1)-y(:,3))) / (constants(2)*(1._pr + constants(1)*sin(y(:,1)-y(:,3))*sin(y(:,1)-y(:,3))))
        
    y(:,4) = ((1._pr + constants(1))*p2 - p1*constants(2)*constants(1)*cos(y(:,1)-y(:,3))) / (constants(2)*constants(2)&
    *constants(1)*(1._pr + constants(1)*sin(y(:,1)-y(:,3))*sin(y(:,1)-y(:,3))))

end subroutine conserve_energy


subroutine get_energy_lagrangian(E, y, constants)
    real(kind=pr), dimension (:), intent(in)                :: constants
    real(kind=pr), dimension(:,:), allocatable, intent(in)  :: y  ! Input: m x 4
    real(kind=pr), dimension(size(y,1),2)                   :: momenta
    real(kind=pr), dimension (:), allocatable, intent(out)  :: E
    
    call get_generalized_momenta(y,constants,momenta)
    E = ((momenta(:,1)*constants(2)-momenta(:,2)*cos(y(:,1)-y(:,3))) **2._pr*2._pr*constants(1)/(2._pr + constants(1) &
    - constants(1)*cos(2._pr*(y(:,1) - y(:,3)))) + momenta(:,2)*momenta(:,2))/(2._pr*constants(2)**2*constants(1))-constants(3)&
    *((1._pr + constants(1))*cos(y(:,1))+constants(2)*constants(1)*cos(y(:,3)))

end subroutine get_energy_lagrangian


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

subroutine generate_flipping_initial_conditions(IC, bounds, N_theta1, N_theta2)
    integer, intent(in)                                     :: N_theta1, N_theta2
    integer                                                 :: ii, jj
    real(kind=pr), dimension(:,:), intent(in)               :: bounds ! 2xm array
    real(kind=pr), dimension(:,:), allocatable, intent(out) :: IC 

    ! This subroutine is grid specific 

    allocate(IC(N_theta1*N_theta2, size(bounds,2)))

    IC = 0
    IC(:,1) = (/((bounds(1,1) + (bounds(2,1) - bounds(1,1))*real(jj,pr)/real(N_theta1-1,pr), jj=0,N_theta1-1), ii=0,N_theta2-1)/)
    IC(:,3) = (/((bounds(1,3) + (bounds(2,3) - bounds(1,3))*real(jj,pr)/real(N_theta2-1,pr), ii=0,N_theta1-1), jj=0,N_theta2-1)/)

end subroutine generate_flipping_initial_conditions

subroutine Hey_chaos(file,f,a,b,h,y0, constants)  
    real (kind=pr), dimension (:,:), allocatable        :: y, y_tmp, y0_tmp, yy0
    real (kind=pr), dimension (:,:), intent (in)        :: y0
    real (kind=pr), dimension (:), intent (in)          :: constants
    logical, dimension (:), allocatable                 :: active_calcs, calcs_tmp
    real (kind=pr), intent (in)                         :: a,b,h
    real (kind=pr)                                      :: t
    integer (kind=int_large)                            :: m, N_inactive
    integer                                             :: ii,jj,kk,NN,unitnum
    character (len=:), allocatable, intent (in)         :: file
    
    interface
    function f(xx,yy, constants)
        real(kind=8), dimension (:,:), allocatable      ::f
        real(kind=8), intent(in)                        :: xx
        real(kind=8), dimension (:), intent(in)         :: constants
        real(kind=8), dimension(:,:), intent(in)        :: yy
    end function
    end interface

    allocate(y(size(y0,1),size(y0,2)), active_calcs(size(y0,1)))
    NN=int((b-a)/h)
    m=size(y0,2)
    N_inactive = 0
    t=a
    y=y0
    yy0 = y0
    active_calcs = .true.


    open(newunit=unitnum,file=file, status = 'replace')

        ! Discard non-flipping systems. Could be generalized by passing another function to indicate the condition for discarding points. Since there are only 2 systems, this simple approach suffices.
        do jj = 1, size(y,1)
            if ((2._pr*cos(y(jj,1)) + cos(y(jj,3))) > 1._pr) then  ! Condition for constants = 1
            !if ((8._pr*cos(y(jj,1)) + cos(y(jj,3))) > 7._pr) then   ! Conditions for alpha=1/3, beta=gamma=1/2
                active_calcs(jj) = .false.
                N_inactive = N_inactive + 1
            end if
        end do
        allocate( y_tmp(size(y,1)-N_inactive,size(y,2)) , y0_tmp(size(yy0,1)-N_inactive,size(yy0,2)))
        allocate(calcs_tmp(size(active_calcs)-N_inactive))
        kk = 1
        do jj = 1, size(y,1)
            if (active_calcs(jj)) then
                y_tmp(kk,:)     =   y(jj,:)
                y0_tmp(kk,:)    =   yy0(jj,:)
                calcs_tmp(kk)   =   active_calcs(jj) 
                kk = kk + 1 
            elseif (.not.(active_calcs(jj))) then
                write(unitnum,format_style) 2._pr*NN*h, yy0(jj,1), yy0(jj,3)
            end if
        end do
        deallocate(y,yy0,active_calcs)
        y  = y_tmp
        yy0 = y0_tmp
        active_calcs = calcs_tmp
        deallocate(y_tmp, y0_tmp, calcs_tmp)
        N_inactive = 0

        ! Calculate the time evolution for the remaining systems and save in file every time more than 100 flips occur
        do ii=1,NN
            call RK_4_array(f,h,y,t,constants)
            t=a+h*ii
            do jj = 1, size(y,1)
                if ((( abs(y(jj,1)) > pi) .or. ( abs(y(jj,3)) > pi)) .and. active_calcs(jj)) then
                    write(unitnum,format_style) t, yy0(jj,1), yy0(jj,3)
                    active_calcs(jj) = .false.
                    N_inactive = N_inactive + 1
                end if
            end do
            if (N_inactive>100) then 
                allocate(y_tmp(size(y,1)-N_inactive,size(y,2)),y0_tmp(size(yy0,1)-N_inactive,size(yy0,2)))
                allocate(calcs_tmp(size(active_calcs)-N_inactive))
                kk = 1
                do jj = 1, size(y,1)
                    if (active_calcs(jj)) then
                        y_tmp(kk,:)     =   y(jj,:)
                        y0_tmp(kk,:)    =   yy0(jj,:)
                        calcs_tmp(kk)   =   active_calcs(jj) 
                        kk = kk + 1 
                    end if
                end do
                deallocate(y,yy0,active_calcs)
                y  = y_tmp
                yy0 = y0_tmp
                active_calcs = calcs_tmp
                deallocate(y_tmp, y0_tmp, calcs_tmp)
                N_inactive = 0
            end if
        end do

        ! For the cases which have not flipped, we save double the time and the initial angles for later plot 
        do jj = 1, size(y,1)
            if (active_calcs(jj)) then
                write(unitnum,format_style) 2._pr*t, yy0(jj,1), yy0(jj,3)
                active_calcs(jj) = .false.
            end if
        end do
    close(unitnum)

    deallocate(y)

end subroutine Hey_chaos

!############################################################################################################## 
! NOT USED: In process of debugging 
!############################################################################################################## 

subroutine Poincare_section_hamiltonian(file,f,a,b,h,y0, constants, energy)   
    real (kind=pr), dimension (:,:), allocatable               :: y, y_tmp
    real (kind=pr), dimension (:,:), intent (in)               :: y0
    real (kind=pr), dimension (:), intent (in)                 :: constants
    real (kind=pr), intent (in)                                :: a,b,h
    real (kind=pr), intent (in)                                :: energy
    real (kind=pr)                                             :: t, angle
    integer (kind=int_large)                                   :: m
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
    m=size(y0,2)
    t=a
    y=y0
    y_tmp=y0

    open(newunit=unitnum,file=file, status = 'replace')
        do ii=1,NN
            call RK_4_array(f,h,y,t,constants)
            t = a + h*ii

    !        do jj = 1, size(y,1)
    !            print*,"P1:",y(jj,2),"P2:",y(jj,4), "Energy:", ((y(jj,2)*constants(2)-y(jj,4)*cos(y(jj,1)-y(jj,3)))&
    !            **2._pr*2._pr*constants(1)/(2._pr + constants(1) - constants(1)*cos(2._pr*(y(jj,1) - y(jj,3)))) + y(jj,4)*y(jj,4))&
    !            /(2._pr*constants(2)**2*constants(1))-constants(3)*((1._pr + constants(1))*cos(y(jj,1))+constants(2)*constants(1)&
    !            *cos(y(jj,3)))
            
    !            if ((2._pr*constants(2)*constants(2)*constants(1)*(energy + constants(3)*((1._pr + constants(1))*cos(y(jj,1))&
    !            +constants(2)*constants(1)*cos(y(jj,3)))) - y(jj,4)*y(jj,4)) < 0._pr) then
    !                print*,"Kinetic energy surpassed the total energy in step:",ii, "for initial condition number:", jj
    !                print*, "Break energy:", ((y(jj,2)*constants(2)-y(jj,4)*cos(y(jj,1)-y(jj,3)))**2._pr*2._pr*constants(1)&
    !                /(2._pr + constants(1) - constants(1)*cos(2._pr*(y(jj,1) - y(jj,3)))) + y(jj,4)*y(jj,4))/(2._pr*constants(2)&
    !                **2*constants(1))-constants(3)*((1._pr + constants(1))*cos(y(jj,1))+constants(2)*constants(1)*cos(y(jj,3)))
                
   !                 print*, "Break conditions:", y(jj,:)
   !             end if
   !         end do

    !        call conserve_energy_hamiltonian(energy, y, constants)

            

            do jj = 1, size(y,1)
                if ((y_tmp(jj,1)*y(jj,1) <= 0.0_pr).and. y(jj,2)>0) then
                    angle = mod(y(jj,3)+pi, 2._pr*pi) - pi
                    angle = mod(angle-pi, 2._pr*pi) + pi
                    write(unitnum,format_style) angle, y(jj,4)
                end if
            end do
            y_tmp = y
        end do
    close(unitnum)

    deallocate(y)

end subroutine Poincare_section_hamiltonian

subroutine conserve_energy_hamiltonian(E, y, constants)
    real(kind=pr), intent(in)                        :: E
    real(kind=pr), dimension (:), intent(in)         :: constants
    real(kind=pr), dimension(:,:), intent(inout)     :: y  ! Input: m x 4
    
    ! Calculate p1 given q1, p2 and q2 such that E is conserved
    y(:,2) = y(:,4)*cos(y(:,1)-y(:,3))/constants(2) + sqrt((2._pr + constants(1) - constants(1)*cos(2._pr*(y(:,1) - y(:,3)))) &
    *(2._pr*constants(2)*constants(2)*constants(1)*(E + constants(3)*((1._pr + constants(1))*cos(y(:,1))+constants(2)*constants(1) &
    *cos(y(:,3)))) - y(:,4)*y(:,4))) /(constants(2)*sqrt(constants(1)*2._pr))

end subroutine conserve_energy_hamiltonian

END MODULE