!***************************************************************
! Program: ex3.f90
! Purpose: simulation of a double pendulum
!
!
! Description: calculates the trajectories of the pendulums (do_trajectories). Calculates the Poincare Calculated the Poincare sections for the conditions such that the angle of the first pendulum is 0 and its momentum is positive (do_poincare_lagrangian and do_poincare_hamiltonian but this last one is not debugged). Calculates the power spectrum with discrete Fourier transform (do_power_spectrum). Calculates the time it takes for a pendulum to flip, given different initial conditions as proposed by Heyl (do_hey). 
!
! Input: dt = time step. constants = constants for the non-dimensional formulation. bounds (:,i) = upper (i=2) or lower (i=1) bounds of the initial conditions, when used, j is used to specify amount of initial conditions to be used. initial_conditions can be used to specify manually the initial conditions.
!
!
! Output: 
!
!
!
! Room for improvement: The subroutine Hey_chaos can be improved by cutting in half the  size of the arrays that save the initial
!   conditions information. In this way the allocation should be even faster and we would use half the memory.  
!                       The system should be adimensionalized so the time evolution is between 0 and 1. To be implemented.
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex3
    use, intrinsic :: iso_c_binding
    use precision
    use funciones
    use subrutinas
    use constantes
    implicit none
    include "/usr/include/fftw3.f03" 

    real(kind=pr)                               :: start_time, end_time, dt, energy, factor
    real(kind=pr), dimension(:,:), allocatable  :: initial_conditions, bounds
    real(kind=pr), dimension(:), allocatable    :: times, constants, angles1, angles2, energies
    real(kind=pr), dimension(:,:), allocatable  :: angular_trajectories, cartesian_trajectories
    integer                                     :: unitnum, initial_condition_index, coordinate
    integer(kind=int_large)                     :: i, j, N_theta1, N_theta2, N_steps, transitory_steps
    character (len=:), allocatable              :: file_angular, file_cartesian, prefix_power_spectrum, file_power_spectrum
    character (len=:), allocatable              :: prefix_poincare , file_poincare, suffix, file_Hey
    logical                                     :: do_trajectories,do_power_spectrum
    logical                                     :: do_poincare_hamiltonian, do_poincare_lagrangian, do_hey
    type(C_ptr)                                 :: plan_rc
    real(C_double), allocatable, dimension(:)   :: sample_in, angles_sample
    complex(C_double_complex), allocatable, dimension(:)    :: fftw_out

    abstract     interface
        function funcion(x_x,y_y, constants)
            use precision
            implicit none
            real(kind=pr), dimension (:,:), allocatable  ::funcion
            real(kind=pr), intent(in)                    :: x_x
            real(kind=pr), dimension (:), intent(in)     :: constants
            real(kind=pr), dimension(:,:), intent(in)    :: y_y
        end function  
    end interface

    procedure (funcion), pointer :: function_pointer => null()
    function_pointer => lagrangian_formalism
        
    start_time      =   0._pr
    end_time        =   10000._pr
    dt              =   0.01_pr
    constants       =   (/1._pr/3._pr,0.5_pr,0.5_pr/)
    file_angular    =   "./datos/datos_ex3_RK4_ang.out"
    file_cartesian  =   "./datos/datos_ex3_RK4_cart.out"
    prefix_poincare =   "./datos/poincare/section_E_"
    suffix          =   ".out"
    file_Hey        =   "./datos/datos_ex3_Hey.out"

    do_trajectories = .False.
    do_power_spectrum  = .False.
    do_poincare_lagrangian  = .False.
    do_hey          = .True.
    do_poincare_hamiltonian = .True. ! Not debugged

!############################################################################################################## 

    if (do_trajectories) then
        allocate(initial_conditions(2,4))
        initial_conditions(1,:)=(/0._pr,0.332_pr,0._pr,0.845_pr/)  !angle1, angle1', angle2, angle2'
        initial_conditions(2,:)=(/0._pr,sqrt(1.125_pr),0._pr,0._pr/)

        ! This commented section was to analyze other energies for the trajectories
        !    initial_conditions (1,:) = (/-0._pr,0._pr,0.0_pr,-0.1_pr/) ! Lower bounds for initial conditions
        !    initial_conditions (2,:) = (/0._pr,0._pr,-0.0_pr,0.1_pr/)  ! Upper bounds for initial conditions
        !    energy = -0.745_pr 
        !    call conserve_energy(energy, initial_conditions, constants)

        initial_condition_index = 1 ! This index specifies that the coordinates associated to said index will be converted to cartesian
        ! Calculate tajectories
        Print*, "--------------------    Calculating trajectories    --------------------------" 
        print*, "Initial conditions: theta1, theta1', theta2, theta2' "
        print*, ""
        call DES_write(file_angular,function_pointer,start_time,end_time,dt,initial_conditions, constants)


        ! Read previous calculations (can be improved) and transform to cartesian coordinates.
        Print*, "-----------    Converting angular coordinates to cartesian    -----------------" 
        print*, ""

        call get_cart_coordinates(file_angular,times,angular_trajectories, cartesian_trajectories, constants&
        , initial_condition_index)
        print*, "Total frames:", size(times)," with dt = ", (times(2) - times(1)), "s and Δt = ", (end_time-start_time)
        
        ! Save to another file. We save x2-x1 and y2-y1 for  later plotting
        open (newunit=unitnum, file=file_cartesian, status='replace')
            write(unitnum,*) "## t | x1 | y1 | x2 | y2 | x2-x1 | y2-y1"
            DO i = 1, size(times)
                write(unitnum, format_style) times(i), cartesian_trajectories(i,1), cartesian_trajectories(i,2), &
                cartesian_trajectories(i,3), cartesian_trajectories(i,4), cartesian_trajectories(i,3)-cartesian_trajectories(i,1), &
                cartesian_trajectories(i,4)-cartesian_trajectories(i,2)   
            END DO
        close(unitnum)
        deallocate(initial_conditions)


    end if

!############################################################################################################## 

    if (do_poincare_lagrangian) then
        function_pointer => lagrangian_formalism

        Print*, "-----------    Calculating Poincare Sections    -----------------" 
        print*, "Initial conditions: theta1, theta1', theta2, theta2'"

        ! Defining lower and upper for initial conditions
        allocate(bounds(2, 4))  
        bounds (1,:) = (/-0.2_pr,0._pr,0.4_pr,-4._pr/) ! Lower bounds for initial conditions
        bounds (2,:) = (/0.2_pr,0._pr,-0.4_pr,4._pr/)  ! Upper bounds for initial conditions
        j=100    ! Number of initial conditions
        call generate_initial_conditions(initial_conditions,bounds,j)

        ! One can loop over the energy values but the bounds could be modified when E changes, thus we only use a scalar energy. Alternatively, we could have a 3D array for the bounds where the 3rd dimension is used to specify the bounds for a given energy in the 1D energy array. Such a complication is unnecessary. 
        energy = -0.0_pr 
        call conserve_energy(energy, initial_conditions, constants)
        do i = 1, size(initial_conditions,1)    
            call get_energy_lagrangian(energies, initial_conditions, constants)
            print*,"Initial Condition:", initial_conditions(i,:)&
            ,"   Energy:",energies(i)
        end do


        call create_file_name(prefix_poincare, energy, suffix, file_poincare)
        call Poincare_section(file_Poincare,function_pointer,start_time,end_time,dt,initial_conditions, constants)
        
        deallocate(initial_conditions,bounds)
    end if

!############################################################################################################## 

    if (do_power_spectrum) then
        transitory_steps = 1        ! Must be chosen greater than 0 
        initial_condition_index = 1 ! This index specified that the coordinates associated to said index will be converted to cartesian
        coordinate = 2              ! Choose which coordinate is to be calculated (theta1 or theta2)
        energy = -0.745_pr          ! See do_poincare_section
        allocate(initial_conditions(2,4))
        initial_conditions(1,:)=(/0._pr,0._pr,0._pr,0.15_pr/)  !angle1, angle1', angle2, angle2'
        initial_conditions(2,:)=(/0._pr,0._pr,0._pr,0.15_pr/)
        call conserve_energy(energy, initial_conditions, constants)
        do i = 1, size(initial_conditions,1)
            print*, "Initail condition N^0", i,":", initial_conditions(i,:)
        end do

        ! Calculate tajectories
        Print*, "--------------------    Calculating trajectories    --------------------------" 
        print*, "Initial conditions: theta1, theta1', theta2, theta2' "
        print*, ""
        call DES_write(file_angular,function_pointer,start_time,end_time,dt,initial_conditions, constants)


        ! Read previous calculations (can be improved) and transform to cartesian coordinates.
        Print*, "-----------    Getting all angular coordinates    -----------------" 
        print*, ""

        call get_ang_coordinates(file_angular,times,angular_trajectories)
        deallocate(initial_conditions)
        print*, "Total steps:", size(times)," with dt = ", (times(2) - times(1)), "s and Δt = ", (end_time-start_time)
        

        Print*, "-----------     Calculating FFT      -----------------" 
        print*, ""

        angles1 = angular_trajectories(:,initial_condition_index)
        angles2 = angular_trajectories(:,initial_condition_index + 4)
        
        select case(coordinate)
            case default
                angles_sample = (/(angles1(i), i = transitory_steps, size(angles1,1))/)
                print*, "Using default coordinate for FFT: theta1"
                prefix_power_spectrum =   "./datos/FFTW/FFTW_theta1_IC_"
            case (1)
                angles_sample = (/(angles1(i), i = transitory_steps, size(angles1,1))/)
                print*, "Using coordinate for FFT: theta1"
                prefix_power_spectrum =   "./datos/FFTW/FFTW_theta1_IC_"
            case (2)
                angles_sample = (/(angles2(i), i = transitory_steps, size(angles1,1))/)
                print*, "Using coordinate for FFT: theta2"
                prefix_power_spectrum =   "./datos/FFTW/FFTW_theta2_IC_"
        end select

        call create_file_name(prefix_power_spectrum, energy, suffix, file_power_spectrum)
        
        
        ! Create FFTW plans
        N_steps = size(angles_sample)
        allocate(sample_in(N_steps), fftw_out(N_steps/2+1))

    
        plan_rc = fftw_plan_dft_r2c_1d(N_steps, sample_in, fftw_out, FFTW_MEASURE)

       
        call fftw_execute_dft_r2c(plan_rc, angles_sample, fftw_out)
        factor = 1._pr/(end_time - start_time)
        call FFTW_forward_write(file_power_spectrum,fftw_out,N_steps, factor)
        call fftw_destroy_plan(plan_rc)

        deallocate(file_power_spectrum, angles_sample, sample_in)
    end if
        
!############################################################################################################## 

    if (do_hey) then 
        Print*, "-----------    Calculating Hey Chaos    -----------------" 
        print*, ""
        !constants(:) =  1._pr
        allocate(bounds(2, 4))  
        bounds (1,:) = (/-3._pr,0._pr,0._pr,0._pr/)
        bounds (2,:) = (/3._pr,0._pr,3._pr,0._pr/)
        N_theta1 = 1200 ! Number of the angles of the first pendulum to be calculated
        N_theta2 = 600  ! Number of the angles of the second pendulum to be calculated

        call generate_flipping_initial_conditions(initial_conditions, bounds, N_theta1, N_theta2)

        call Hey_chaos(file_Hey,function_pointer,start_time,end_time,dt,initial_conditions, constants)
    end if

!############################################################################################################## 
! NOT USED:   THIS SECTION IS DONE IN THE HAMILTONIAN FORMALISM!!!
!############################################################################################################## 
    
    if (do_poincare_hamiltonian) then
        function_pointer => hamiltonian_formalism

        Print*, "-----------    Calculating Poincare Sections    -----------------" 
        print*, "Initial conditions: theta1, P1, theta2, P2"

        ! Defining lower and upper for initial conditions
        allocate(bounds(2, 4))  
        bounds (1,:) = (/0._pr,0._pr,0._pr,-0.35_pr/) ! Lower bounds for initial conditions
        bounds (2,:) = (/0._pr,0._pr,0._pr,0.35_pr/)  ! Upper bounds for initial conditions
        j=50    ! Number of initial conditions
        call generate_initial_conditions(initial_conditions,bounds,j)

        !One can loop over the energy values but the bounds could be modified when E changes
        energy = -0.0_pr 
        call conserve_energy_hamiltonian(energy, initial_conditions, constants)
        do i = 1, size(initial_conditions,1) 
            print*, "Initial Condition:", initial_conditions(i,:)&
            ,"   Energy:", ((initial_conditions(i,2)*constants(2)-initial_conditions(i,4)*cos(initial_conditions(i,1)&
            -initial_conditions(i,3)))**2._pr*2._pr*constants(1)/(2._pr + constants(1) - constants(1)*cos(2._pr&
            *(initial_conditions(i,1) - initial_conditions(i,3)))) + initial_conditions(i,4)*initial_conditions(i,4))&
            /(2._pr*constants(2)**2*constants(1))-constants(3)*((1._pr + constants(1))*cos(initial_conditions(i,1))+constants(2)&
            *constants(1)*cos(initial_conditions(i,3)))
        end do


        call create_file_name(prefix_poincare, energy, suffix, file_poincare)
        call Poincare_section_hamiltonian(file_Poincare,function_pointer,start_time,end_time,dt,initial_conditions&
        , constants, energy)
        
        deallocate(initial_conditions,bounds)
    end if


end program ex3
