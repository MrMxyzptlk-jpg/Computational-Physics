MODULE subrutinas
    use, intrinsic :: iso_c_binding
    use precision
    use funciones
    implicit none

    CONTAINS

subroutine FFTW_forward_write(factor,unit_pow,file_pow,fftw_out,N)
    integer                                     :: N, i, unit_pow
    real(pr)                                    :: factor 
    character(len=:), allocatable               :: file_pow
    complex(C_double_complex), allocatable, dimension(:)    :: fftw_out
    
    open(unit=unit_pow, file=file_pow, status='replace')
        write(unit_pow, *) "# Frequency | FFTW (Re) | FFTW (Im)"
        do i = -N/2, N/2
            if (i < 0) then
                write(unit_pow, *) factor * real(i, pr), real(conjg(fftw_out(-i+1)) / real(N, pr)), &
                aimag(conjg(fftw_out(-i+1)) / real(N, pr))
            else
                write(unit_pow, *) factor * real(i, pr), real(fftw_out(i+1) / real(N, pr)), aimag(fftw_out(i+1) / real(N, pr))
            end if
        end do
        print*,  '            Values saved in: ', file_pow
        print*, ''
    close(unit_pow)
END SUBROUTINE FFTW_forward_write

    
subroutine FFTW_inverse_write(step,unit_inv,file_inv, sample_in, sample_out,N)
    integer                                     :: N, i, unit_inv
    real(pr)                                    :: step 
    character(len=:), allocatable               :: file_inv
    real(C_double), allocatable, dimension(:)   :: sample_in, sample_out
    
    open(unit=unit_inv, file=file_inv, status='replace')
        write(unit_inv, *) "# Position | Sampled Function | FFTW Inverse Transform "
        do i = 1, N
            write(unit_inv, *) step*(i-1), sample_in(i), sample_out(i) / real(N, pr)
        end do
        print*,  '            Values saved in: ', file_inv
        print*, ''
    close(unit_inv)
END SUBROUTINE FFTW_inverse_write

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

subroutine logistic_map_step(position, parameter_r, index)
    real(kind=pr),intent(in),dimension(:)   :: parameter_r
    real(kind=pr),intent(inout),dimension(:):: position
    integer,intent(in)                      :: index
    
    position = parameter_r(index)*position*(1._pr-position)
END SUBROUTINE logistic_map_step

subroutine logistic_map(N_steps,initial_conditions, parameter_r, prefix, suffix)
    integer                                         :: unitnum, step, index
    integer, intent(in)                             :: N_steps
    real(pr), allocatable, dimension(:), intent(in) :: initial_conditions, parameter_r
    real(pr), allocatable, dimension(:)             :: position
    character(len=:), allocatable                   :: filename, suffix, prefix

    Print*, "-----------     Calculating trajectories of the logistic map     -----------------" 
    print*, ""

    ! Loop through all parameters and calculate all trayectories simultaneously
    DO index = 1, size(parameter_r)
        ! Calculating the first 300 irrelevant steps
        position = initial_conditions
        DO step = 1, 300
            call logistic_map_step(position, parameter_r, index)
        END DO

        ! Create the ./datos/trayectory_IC_*****_param_*****.out file name
        call create_file_name(prefix,parameter_r(index),suffix,filename) 

        ! Calculate and save the trayectory
        !print*, "Calculating trajectories with ", N_steps, " steps and r = ", parameter_r(index)
        open(unit=unitnum, file=filename, status='replace')
            !write(unitnum,*)"# Step | Positions | Parameter "
            DO step = 1, N_steps
                call logistic_map_step(position, parameter_r, index)
                write(unitnum,*) step, position , parameter_r(index)
            END DO
        close(unitnum)
        !print*, "Steps and positions written to file: ", filename
        deallocate(filename)
    END DO
    print*, ""


END subroutine logistic_map


subroutine logistic_map_and_lyapunov(N_steps,initial_conditions, parameter_r, prefix, suffix, file_lyapunov)
    integer                                         :: unitnum,unitnum2, step, index
    integer, intent(in)                             :: N_steps
    real(pr)                                        :: lyapunov_scalar
    real(pr), allocatable, dimension(:), intent(in) :: initial_conditions, parameter_r
    real(pr), dimension(size(initial_conditions))   :: position, lyapunov_exponent
    character(len=:), allocatable                   :: filename, suffix, prefix, file_lyapunov


    Print*, "-----------     Calculating trajectories of the logistic map and Lyapunov exponent     -----------------" 
    print*, ""

    ! Loop through all parameters and calculate all trayectories simultaneously
        open(newunit=unitnum2, file=file_lyapunov, status='replace')
    DO index = 1, size(parameter_r)
        lyapunov_exponent = 0._pr
        ! Calculating the first 300 irrelevant steps
        position = initial_conditions
        DO step = 1, 300
            call logistic_map_step(position, parameter_r, index)
        END DO

        ! Create the ./datos/trayectory_IC_*****_param_*****.out file name
        call create_file_name(prefix,parameter_r(index),suffix,filename) 

        ! Calculate and save the trayectory
        !print*, "Calculating trajectories with ", N_steps, " steps and r = ", parameter_r(index)
        open(newunit=unitnum, file=filename, status='replace')
            !write(unitnum,*)"# Step | Positions | Parameter "
            DO step = 1, N_steps
                call logistic_map_step(position, parameter_r, index)
                lyapunov_exponent = lyapunov_exponent + log(abs(parameter_r(index)-2._pr*parameter_r(index)*position))
                write(unitnum,*) step, position , parameter_r(index)
            END DO
        close(unitnum)
        deallocate(filename)

        lyapunov_scalar = sum(lyapunov_exponent)/real(size(lyapunov_exponent)*N_steps,pr)
            write(unitnum2,*) index, lyapunov_scalar, parameter_r(index)
    END DO
    close(unitnum2)
    print*, ""

END subroutine logistic_map_and_lyapunov


END MODULE
