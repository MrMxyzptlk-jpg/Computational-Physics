!***************************************************************
! Program: ex2.f90
! Purpose: 
!
!
! Description: Calculate the trajectories for different cases of the logistic equation.
!   Calculate the FFT for a single initial conditions but different values of the logistic parameter.
!   Create a normalized histogram for the trajectory starting from a single initial condition with 
!   different values of the logistic parameter. Create 'orbit diagrams' (i.e. bifurcation diagrams).
!   Calculate and graph the Lyapunov exponent.
!
! Input: Number of steps, set of initial conditions, set of logistic parameters.
!
! Output: Set of files in their corresponding folders inside the './datos' folder, ready to graph using the 
!   plot.gp gnuplot script. 
!
!
! Room for improvement: The main subroutine utilized is the 'logistic map' which uses arrays for 
!   different initial conditions but does loops for different values of the logistic parameter. This 
!   could be improved by using nxm dimensional arrays with n the number of initial conditions and m the
!   number of different parameters. That being said, for this small calculations the loop works well enough.
!
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex2
    use, intrinsic :: iso_c_binding
    use precision
    use funciones
    use subrutinas
    implicit none
    include "/usr/include/fftw3.f03" 


    integer                                     :: j, unitnum, N_steps, step, index
    real(pr), allocatable,dimension(:)          :: initial_conditions, parameter_r, trajectory, array_tmp
    character(len=:), allocatable               :: filename, suffix, prefix, file_pow,file_lyapunov
    real(pr)                                    :: lower_limit, upper_limit, range, factor
    logical                                     :: ex2_a, ex2_b, ex2_c, ex2_d, ex2_d_extra, ex2_e
    type(C_ptr)                                 :: plan_rc
    real(C_double), allocatable, dimension(:)   :: sample_in
    complex(C_double_complex), allocatable, dimension(:)    :: fftw_out


    ! Define all relevant information
    N_steps = 1000       ! Number of relevant steps (300 extra steps will be done first, and discarded)
    initial_conditions = [0.06_pr, 0.03_pr, 0.6_pr, 0.9_pr]
    parameter_r = [1.5_pr, 3.3_pr, 3.5_pr, 3.55_pr, 4._pr]

    ex2_a = .False.
    ex2_b = .False.
    ex2_c = .False.
    ex2_d = .False.
    ex2_d_extra = .False.
    ex2_e = .False.
    
!*****************************************************************************************************
!       END OF SETTINGS SECTION / START OF THE PROGRAM
!*****************************************************************************************************

    prefix = "./datos/trajectories/trajectories_param_"
    suffix = ".out"
    file_pow = "./datos/FFTW/FFTW_power_spectrum_ex2_param_"

    
    !***********************    Trajectories (ex2-a)   **********************   

    if(ex2_a.or.ex2_b)then
        ! Loop through all parameters and calculate all trayectories simultaneously
        call  logistic_map(N_steps,initial_conditions, parameter_r, prefix, suffix)
    end if

    !***********************    FFTW  (ex2-b) **********************   
    if(ex2_b)then
        ! Create FFTW plans
        allocate(sample_in(N_steps), fftw_out(N_steps/2+1))

        allocate(trajectory(N_steps), array_tmp(size(initial_conditions)))

        ! Read previous calculations with the initial condition x=0.6 and calculate the FFT of each trajectory
        Print*, "-----------     Calculating FFT for the initial condition x=0.6     -----------------" 
        print*, ""
        DO index = 1, size(parameter_r)
            plan_rc = fftw_plan_dft_r2c_1d(N_steps, sample_in, fftw_out, FFTW_MEASURE)
            
            call create_file_name(prefix,parameter_r(index),suffix,filename) 
            open (unit=unitnum, file=filename, status='old',    &
                access='sequential', form='formatted', action='read' )
                DO j = 1, N_steps
                    read (unitnum, *)  step, array_tmp   !positions given different initial conditions (ic)
                    trajectory(j) = array_tmp(3)
                END DO
            close(unitnum)
            deallocate(filename)


            call create_file_name(file_pow,parameter_r(index),suffix,filename) 
            
            call fftw_execute_dft_r2c(plan_rc, trajectory, fftw_out)
            factor = 1._pr / real(N_steps,pr)
            call FFTW_forward_write(factor,unitnum,filename,fftw_out,N_steps)
            call fftw_destroy_plan(plan_rc)
            deallocate(filename)
        END DO

        deallocate(trajectory, array_tmp)
    end if


    !***********************    Trajectories (ex2-c)   **********************    
    if(ex2_c)then
        ! Define all relevant information
        N_steps = 10000       ! Number of relevant steps (300 extra steps will be done first, and discarded)
        initial_conditions = [0.6_pr]
        parameter_r = [4._pr]
        prefix = "./datos/histogram_data_param_"
        suffix = ".out"
        
        call  logistic_map(N_steps,initial_conditions, parameter_r, prefix, suffix)

        deallocate(parameter_r,initial_conditions, prefix)
    end if

    !***********************    Trajectories (ex2-d)   **********************    
    if(ex2_d)then
        ! Define all relevant information
        N_steps = 4096      ! Number of relevant steps (300 extra steps will be done first, and discarded)
        initial_conditions = [0.6_pr]
        prefix = "./datos/orbits/orbit_param_"
        suffix = ".out"

        ! Defining range for r
        lower_limit = 3.4_pr
        upper_limit = 4.0_pr
        range = (upper_limit - lower_limit )
        step = 300
        parameter_r = (/(lower_limit + real(j-1, pr)*range/real(step,pr),j=1,step)/)

        ! Too many files written... could be improved as mentioned at the start of this file 
        call  logistic_map(N_steps,initial_conditions, parameter_r, prefix, suffix)

        if(ex2_d_extra)then ! This should have been implemented in a submodule for tidiness (it's the same as before)
            ! Create FFTW plans
            allocate(sample_in(N_steps), fftw_out(N_steps/2+1))
    
            allocate(trajectory(N_steps), array_tmp(size(initial_conditions)))
    
            ! Read previous calculations with the initial condition x=0.6 and calculate the FFT of each trajectory
            Print*, "-----------     Calculating FFT for the initial condition x=0.6     -----------------" 
            print*, ""
            DO index = 1, size(parameter_r)
                plan_rc = fftw_plan_dft_r2c_1d(N_steps, sample_in, fftw_out, FFTW_MEASURE)
                call create_file_name(prefix,parameter_r(index),suffix,filename) 
                open (newunit=unitnum, file=filename, status='old',    &
                    access='sequential', form='formatted', action='read' )
                    DO j = 1, N_steps
                        read (unitnum, *)  step, array_tmp   !positions given different initial conditions (ic)
                        trajectory(j) = array_tmp(1)
                    END DO
                close(unitnum)
                deallocate(filename)
                
                factor = 1._pr / real(N_steps,pr)
    
                call create_file_name(file_pow,parameter_r(index),suffix,filename) 
                
                call fftw_execute_dft_r2c(plan_rc, trajectory, fftw_out)
                call FFTW_forward_write(factor,unitnum,filename,fftw_out,N_steps)
                call fftw_destroy_plan(plan_rc)
            END DO
            deallocate(trajectory, array_tmp)
            print*, ""
        end if


        deallocate(parameter_r)

        N_steps = 1200       ! Number of relevant steps (300 extra steps will be done first, and discarded)
        
        ! Defining range for r
        lower_limit = 3.847_pr
        upper_limit = 3.8568_pr
        range = (upper_limit - lower_limit )
        step = 200
        parameter_r = (/(lower_limit + real(j-1, pr)*range/real(step,pr),j=1,step)/)

        ! Too many files written... could be improved as mentioned at the start or this file 
        call  logistic_map(N_steps,initial_conditions, parameter_r, prefix, suffix)
        
        deallocate(parameter_r,initial_conditions, prefix)
    end if


 !***********************    Trajectories (ex2-e)   **********************    
    if(ex2_e)then
        ! Define all relevant information
        N_steps = 500      ! Number of relevant steps (300 extra steps will be done first, and discarded)
        initial_conditions = [0.6_pr]
        prefix = "./datos/lyapunov/orbit_param_"
        suffix = ".out"
        file_lyapunov = "./datos/lyapunov/lyapunov_exponents.out"

        ! Defining range for r
        lower_limit = 2.001_pr
        upper_limit = 4.0_pr
        range = (upper_limit - lower_limit )
        step = 500
        parameter_r = (/(lower_limit + real(j-1, pr)*range/real(step,pr),j=1,step)/)

        ! Too many files written... could be improved as mentioned at the start of this file 
        call  logistic_map_and_lyapunov(N_steps,initial_conditions, parameter_r, prefix, suffix, file_lyapunov)

        deallocate(parameter_r,initial_conditions, prefix)
    end if


end program ex2
