!***************************************************************
! Program: ex4.f90
! Purpose: analisis of the Pullen-Edmonds Hamiltonian
!
!
! Description:  Analises the convergence of the time the time step and the drift in energy due to the time evolution (do_convergence)calculate the dependence of the coordinates and momentums with time (do_trajectories). Calculates the power spectrum with discrete Fourier for each coordinate and momentum (do_power_spectrum). Calculated the Poincare sections for the conditions such that one coordinate is 0 and its momentum is positive.
!
! Input:
!
!
! Output: 
!
!
!
! Room for improvement: 
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex4
    use, intrinsic :: iso_c_binding
    use precision
    use funciones
    use subrutinas
    use constantes
    implicit none
    include "/usr/include/fftw3.f03" 

    real(kind=pr)                               :: start_time, end_time, h, energy, factor
    real(kind=pr), dimension(:,:), allocatable  :: initial_conditions, bounds, y, trajectories
    real(kind=pr), dimension(:), allocatable    :: constants, energies, times, errors
    character (len=:), allocatable              :: prefix_trajectory
    character (len=:), allocatable              :: prefix_poincare, prefix_power_spectrum, suffix, filename, file_power_spectrum 
    logical                                     :: do_trajectories, do_power_spectrum, do_poincare 
    logical                                     :: do_convergence
    integer                                     :: i, j, unitnum, io, N_steps, lines, coordinate
    real(kind=pr), dimension (:), allocatable   :: energy_tmp, energy_0_tmp, energy_error
    type(C_ptr)                                 :: plan_rc
    real(C_double), allocatable, dimension(:)   :: sample_in
    complex(C_double_complex), allocatable, dimension(:)    :: fftw_out


    abstract interface
        function funcion(x_x,y_y, constants)
            use precision
            implicit none
            real(kind=pr), dimension (:,:), allocatable             :: funcion
            real(kind=pr), intent(in)                               :: x_x
            real(kind=pr), dimension (:), intent(in)                :: constants
            real(kind=pr), dimension(:,:), intent(in)               :: y_y
        end function
    end interface

    procedure (funcion), pointer :: function_pointer => null()
    function_pointer => hamiltonian_formalism
        
    start_time      =   0._pr
    end_time        =   200._pr
    h               =   0.001_pr
    constants       =   (/0.05_pr/)
    prefix_trajectory  =   "./datos/trajectories/E_"
    prefix_poincare =   "./datos/poincare/section_E_"
    prefix_power_spectrum =   "./datos/FFTW/FFTW_IC_"
    suffix          =   ".out"


    do_convergence      = .False.
    do_trajectories     = .False.
    do_poincare         = .False.
    do_power_spectrum   = .True.

    
!############################################################################################################## 

    if(do_convergence) then
        Print*, "------------------    Calculating trajectories    ------------------------" 
        print*, ""

        ! Defining lower and upper for initial conditions
        allocate(bounds(2, 4))  
        bounds (1,:) = (/-0._pr,0._pr,-1._pr,-3.16_pr/)
        bounds (2,:) = (/0._pr,0._pr,1._pr,3.16_pr/)
        j=1
        call generate_initial_conditions(initial_conditions,bounds,j)

        !One can loop over the energy values but the bounds could be modified when E changes
        energy = 100._pr
        call conserve_energy(energy, initial_conditions, constants) 
        do i = 1, size(initial_conditions,1)
            print*, initial_conditions(i,:)
        end do


        Print*, "--------------------    Calculating trajectories    --------------------------" 
        print*, ""
        !call DES_write(filename,function_pointer,start_time,end_time,h,initial_conditions, constants)

        filename = "./datos/dt_convergence.out"
        allocate(energy_0_tmp(1))
        energy_0_tmp = energy
        open(newunit=unitnum,file=filename,status='replace')
            write(unitnum,*) "## dt | E"
            do i = 1, 10000
                h = 1e-5*i
                call DES(function_pointer,start_time,end_time,h,initial_conditions, y, constants)
                call get_energy(energy_tmp,y, constants)
                call get_energy_error(energy_0_tmp, energy_tmp, energy_error)
                write(unitnum, *) i, h, energy_error
            end do
        close(unitnum)
        deallocate(initial_conditions, filename)
    end if

!############################################################################################################## 

    if(do_trajectories)then
        allocate(initial_conditions(1, 4))
        Print*, "--------------------    Calculating trajectories    --------------------------" 
        print*, ""
        energies = (/5._pr,20._pr,100._pr/)
        do i = 1, size(energies)
            energy = energies(i)
            initial_conditions(1,:) = (/(2._pr), (0._pr), (0._pr), (sqrt(2._pr*energy-4._pr))/)
            call create_file_name(prefix_trajectory, energy, suffix, filename)
            call DES_write(filename,function_pointer,start_time,end_time,h,initial_conditions, constants)
            deallocate(filename)
        end do 
        deallocate(initial_conditions, energies)
    end if

!############################################################################################################## 

    if(do_poincare)then
        Print*, "------------------    Calculating Poincaré Section    ------------------------" 
        print*, ""

        ! Defining lower and upper for initial conditions
        allocate(bounds(2, 4))  
        bounds (1,:) = (/-0._pr,0._pr,-1._pr,-14._pr/)
        bounds (2,:) = (/0._pr,0._pr,1._pr,14._pr/)
        j=100
        call generate_initial_conditions(initial_conditions,bounds,j)

        !One can loop over the energy values but the bounds could be modified when E changes
        energy = 100._pr
        call conserve_energy(energy, initial_conditions, constants) 
        do i = 1, size(initial_conditions,1)
            print*, initial_conditions(i,:)
        end do


        call create_file_name(prefix_poincare, energy, suffix, filename)
        call Poincare_section(filename,function_pointer,start_time,end_time,h,initial_conditions, constants)
        
        deallocate(initial_conditions, filename)
    end if

!############################################################################################################## 
    
    if (do_power_spectrum) then
        energy = 100._pr
        Print*, "-----------     Reading trajectories for orbit with energy E = ",energy,"     -----------------" 
        print*, ""

        call create_file_name(prefix_trajectory, energy, suffix, filename)
        open (newunit=unitnum, file=filename, status='old', access='sequential', form='formatted', action='read' )
            lines = 0 
            DO
                read(unitnum,*,iostat=io)
                if (io/=0) EXIT
                lines = lines + 1
            END DO
            allocate(trajectories(lines,4), times(lines))
            
            rewind(unitnum)
            DO i =1,lines
                read(unitnum,format_style) times(i),trajectories(i,:), errors
            END DO
        close(unitnum)

        

        Print*, "-----------     Calculating FFT      -----------------" 
        print*, ""
        
        ! Create FFTW plans
        N_steps = size(trajectories,1)
        allocate(sample_in(N_steps), fftw_out(N_steps/2+1))

        factor = 1._pr/(end_time - start_time)
        plan_rc = fftw_plan_dft_r2c_1d(N_steps, sample_in, fftw_out, FFTW_MEASURE)
    
        do coordinate = 1, 4
            select case(coordinate)
                case default
                    sample_in = trajectories(:,1)
                    print*, "Using default coordinate for FFT: q1"
                    prefix_power_spectrum =   "./datos/FFTW/FFTW_q1_E_"
                case (1)
                    sample_in = trajectories(:,1)
                    print*, "Using coordinate for FFT: q1"
                    prefix_power_spectrum =   "./datos/FFTW/FFTW_q1_E_"
                case (2)
                    sample_in = trajectories(:,2)
                    print*, "Using coordinate for FFT: p1"
                    prefix_power_spectrum =   "./datos/FFTW/FFTW_p1_E_"
                case (3)
                    sample_in = trajectories(:,3)
                    print*, "Using coordinate for FFT: q2"
                    prefix_power_spectrum =   "./datos/FFTW/FFTW_q2_E_"
                case (4)
                    sample_in = trajectories(:,4)
                    print*, "Using coordinate for FFT: p2"
                    prefix_power_spectrum =   "./datos/FFTW/FFTW_p2_E_"
            end select

            call create_file_name(prefix_power_spectrum, energy, suffix, file_power_spectrum)
            call fftw_execute_dft_r2c(plan_rc, sample_in, fftw_out)
            call FFTW_forward_write(file_power_spectrum,fftw_out,N_steps, factor)

        end do

        call fftw_destroy_plan(plan_rc)

        deallocate(file_power_spectrum, filename, sample_in)
    end if
        

end program ex4
