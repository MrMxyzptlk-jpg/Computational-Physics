MODULE subrutinas
use precision
use constantes
use formats
use funciones
implicit none

contains


subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                        :: fmt  ! Format descriptor
    character(len=12)                       :: x1   ! Temporary string for formatted real
    character(len=:), allocatable           :: prefix, suffix, filename
    integer(kind=int_huge), intent(in)      :: num  ! Input real number

    fmt = '(I10)'  ! Format integer
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

subroutine tridiagonal_matrix_solver(X, diag_low, diag, diag_up)
    real (kind=pr), dimension (:), intent (inout)       :: X
    real (kind=pr), dimension (size(X)-1), intent (in)  :: diag_low, diag_up
    real (kind=pr), dimension (size(X)), intent (in)    :: diag
    real (kind=pr), dimension (size(diag_up))           :: h
    real (kind=pr), dimension (size(X))                 :: p
    integer (kind=int_large)                            :: i

    h(1) = diag_up(1)/diag(1)
    p(1) = X(1)/(diag(1))

    do i = 2, size(X)-1
        if (abs(diag(i) - diag_low(i-1)*h(i-1)) < tiny(1._pr)) then
            print *, "ERROR in Thomas solver: near-zero pivot at i=", i
            stop
        end if
        h(i) =  diag_up(i)/(diag(i) - diag_low(i-1)*h(i-1))
        p(i) =  (X(i) - diag_low(i-1)*p(i-1))/(diag(i) - diag_low(i-1)*h(i-1))
    end do

    p(size(p)) =  (X(size(p)) - diag_low(size(p)-1)*p(size(p)-1))/(diag(size(p)) - diag_low(size(p)-1)*h(size(p)-1))

    X(size(X)) = p(size(X))

    do i = size(X)-1, 1, -1
        X(i) = p(i) - h(i)*X(i+1)
    end do

end subroutine tridiagonal_matrix_solver

subroutine write_to_file(Y, time, dx, conversion_factors, unitnum)
    real (kind=pr), dimension (:), intent (in)  :: Y, conversion_factors
    real (kind=pr), intent (in)                 :: dx, time
    integer (kind=int_medium)                   :: unitnum
    integer                                     :: i

    do i = 1, size(Y)
        write(unitnum, format_style2) time*conversion_factors(2) &
        , dx*real(i-1,pr)*conversion_factors(1), Y(i)
    end do

end subroutine write_to_file

subroutine check_stability_criteria(N, dt, stability)
    integer             :: N
    real (kind=pr)      :: dt, weight
    logical             :: stability

    weight = dt*(real(N-1,pr))*(real(N-1,pr))  ! weight = dt/(dx*dx)  We reformulate it to avoid division

    if (weight==0.5)then
        print*, "Stability criterion is in the edge: dt/dx^2 = ", weight, "= 1/2"
        stability = .true.
    else if(weight<0.5) then
        print*, "Stability criterion: dt/dx^2 = ", weight, "< 1/2"
        stability = .true.
    else
        print*, "Stability criteria is not satisfied: dt/dx^2 = ", weight, "> 1/2"
        stability = .false.
    end if

end subroutine check_stability_criteria

subroutine forward_euler_Neumann_step(Y, dt)
    real (kind=pr), dimension (:), allocatable, intent (inout)  :: Y
    real (kind=pr), dimension (size(Y))                         :: Y_tmp
    real (kind=pr), intent (in)                                 :: dt
    real (kind=pr)                                              :: weight
    integer (kind=int_large)                                    :: j

    weight = dt*real(size(Y)-1,pr)*real(size(Y)-1,pr)  ! weight = dt/(dx*dx)  We reformulate it to avoid division

    Y_tmp(1) = Y(1)*(1._pr - 2._pr*weight) + 2._pr*weight*Y(2)
    Y_tmp(2:size(y_tmp)-1) = (/(Y(j-1)*weight + Y(j)*(1._pr - 2._pr*weight) + Y(j+1)*weight, j=2,size(Y)-1)/)
    Y_tmp(size(Y_tmp)) = 2._pr*weight*Y(size(Y)-1) + Y(size(Y))*(1._pr-2._pr*weight)
    Y = Y_tmp

end subroutine forward_euler_Neumann_step

subroutine forward_euler_Neumann(initial_conditions, Y, dt, t_min, t_max, filename, conversion_factors,epochs)
    real (kind=pr), dimension (:), intent (in)                  :: initial_conditions
    real (kind=pr), dimension (:), allocatable, intent (out)    :: Y
    real (kind=pr), dimension (:), allocatable                  :: conversion_factors
    real (kind=pr), intent (in)                                 :: dt, t_max,t_min
    real (kind=pr)                                              :: dx
    integer (kind=int_large)                                    :: j
    integer (kind=int_medium), intent(in)                       :: epochs
    integer (kind=int_medium)                                   :: unitnum
    character (len=:), allocatable, intent(in)                  :: filename
    logical                                                     :: write_last_iter

    dx = 1._pr/real(size(initial_conditions)-1,pr)   ! dx = 1/(grid points -1)

    Y =  initial_conditions

    ! Check if all iterations will be written or the last one needs to be saved separately
    if(mod(int((t_max - t_min)/dt),epochs)==0) then
        write_last_iter = .false.
    else
        write_last_iter = .true.
    end if

    open(newunit=unitnum, file=filename, status='replace')
        write(unitnum, *) "# dt [s] = ", dt*conversion_factors(2), "dx [m] = ", dx*conversion_factors(1)
    write(unitnum, *) "# time | x | Temperature "
        call write_to_file(Y, t_min, dx, conversion_factors, unitnum)
        do j = 1, int((t_max - t_min)/dt)
            call forward_euler_Neumann_step(Y, dt)
            ! Save state every 'epochs' iterations
            if (mod(j,epochs) == 0) then
                call write_to_file(Y, t_min + dt*real((j), pr), dx, conversion_factors, unitnum)
            end if
        end do

        ! Save the last iteration if it has not been already saved
        if (write_last_iter) then
            call write_to_file(Y, t_min + dt*real((int((t_max - t_min)/dt)), pr), dx, conversion_factors, unitnum)
        end if
    close(unitnum)

end subroutine forward_euler_Neumann

subroutine backward_euler_Neumann(initial_conditions, Y, dt, t_min, t_max, filename, conversion_factors,epochs)
    real (kind=pr), dimension (:), intent (in)                  :: initial_conditions
    real (kind=pr), dimension (:), allocatable, intent (out)    :: Y
    real (kind=pr), dimension (size(initial_conditions))        :: diag
    real (kind=pr), dimension (size(initial_conditions)-1)      :: diag_low, diag_up
    real (kind=pr), dimension (:), intent(in)                   :: conversion_factors
    real (kind=pr), intent (in)                                 :: dt, t_max,t_min
    real (kind=pr)                                              :: dx, weight
    integer (kind=int_large)                                    :: j
    integer (kind=int_medium), intent(in)                       :: epochs
    integer (kind=int_medium)                                   :: unitnum
    character (len=:), allocatable, intent(in)                  :: filename
    logical                                                     :: write_last_iter

    dx = 1._pr/real(size(initial_conditions)-1,pr)   ! dx = 1/(grid points -1)
    weight = dt*real(size(initial_conditions)-1,pr)*real(size(initial_conditions)-1,pr)  ! weight = dt/(dx*dx)  We reformulate it to avoid division

    Y = initial_conditions
    diag_up(1) = -2._pr*weight
    diag_up(2:size(diag_up))    = -weight
    diag = 1._pr + 2._pr*weight
    diag_low(1:size(diag_low)-1) = -weight
    diag_low(size(diag_low))     = -2._pr*weight

    ! Check if all iterations will be written or the last one needs to be saved separately
    if(mod(int((t_max - t_min)/dt),epochs)==0) then
        write_last_iter = .false.
    else
        write_last_iter = .true.
    end if

    open(newunit=unitnum, file=filename, status='replace')
        write(unitnum, *) "# dt [s] = ", dt*conversion_factors(2), "dx [m] = ", dx*conversion_factors(1)
        write(unitnum, *) "# time | x | Temperature "
        call write_to_file(Y, t_min, dx, conversion_factors, unitnum)

        do j = 1, int((t_max - t_min)/dt)
            call tridiagonal_matrix_solver(Y, diag_low, diag, diag_up)

            ! Save state every 'epochs' iterations
            if (mod(j,epochs) == 0) then
                call write_to_file(Y, t_min + dt*real((j), pr), dx, conversion_factors, unitnum)
            end if
        end do

        ! Save the last iteration if it has not been already saved
        if (write_last_iter) then
            call write_to_file(Y, t_min + dt*real((int((t_max - t_min)/dt)), pr), dx, conversion_factors, unitnum)
        end if
    close(unitnum)

end subroutine backward_euler_Neumann

subroutine crank_nicolson_Neumann(initial_conditions, Y, dt, t_min, t_max, filename, conversion_factors,epochs)
    real (kind=pr), dimension (:), intent (in)                  :: initial_conditions
    real (kind=pr), dimension (size(initial_conditions)),  intent (out)    :: Y
    real (kind=pr), intent (in)                                 :: dt, t_max,t_min
    character (len=:), allocatable, intent(in)                  :: filename
    real (kind=pr), dimension (size(initial_conditions))        :: Y_tmp, diag
    real (kind=pr), dimension (size(initial_conditions)-1)      :: diag_low, diag_up
    real (kind=pr), dimension (:), intent(in)                   :: conversion_factors
    real (kind=pr)                                              :: dx, weight, weight_inverse, diag_element
    integer (kind=int_medium), intent(in)                       :: epochs
    integer (kind=int_medium)                                   :: unitnum
    integer (kind=int_large)                                    :: i, j
    logical                                                     :: write_last_iter

    dx = 1._pr/real(size(initial_conditions)-1,pr)   ! dx = 1/(grid points -1)
    weight = dt*real(size(Y)-1,pr)*real(size(Y)-1,pr)  ! weight = dt/(dx*dx)  We reformulate it to avoid division
    weight_inverse = 1._pr/weight

    Y = initial_conditions
    diag_up(1) = -2._pr
    diag_up(2:size(diag_up))    = -1._pr
    diag = 2._pr + 2._pr*weight_inverse
    diag_low(1:size(diag_low)-1) = -1._pr
    diag_low(size(diag_low))     = -2._pr

    ! Check if all iterations will be written or the last one needs to be saved separately
    if(mod(int((t_max - t_min)/dt),epochs)==0) then
        write_last_iter = .false.
    else
        write_last_iter = .true.
    end if

    diag_element = 2._pr*weight_inverse - 2._pr     ! For the calculation of a temporary matrix in order to solve a tri-diagonal system of equations

    open(newunit=unitnum, file=filename, status='replace')
        write(unitnum, *) "# dt [s] = ", dt*conversion_factors(2), "dx [m] = ", dx*conversion_factors(1)
        write(unitnum, *) "# time | x | Temperature "
        call write_to_file(Y, t_min, dx, conversion_factors, unitnum)
        do j = 1, int((t_max - t_min)/dt)
            Y_tmp(1) = Y(1)*diag_element + Y(2)*2._pr
            Y_tmp(2:size(y_tmp)-1) = (/(Y(i-1) + Y(i)*diag_element + Y(i+1), i=2, size(Y)-1)/)
            Y_tmp(size(y_tmp)) = 2._pr*Y(size(Y)-1) + Y(size(Y))*diag_element

            call tridiagonal_matrix_solver(Y_tmp, diag_low, diag, diag_up)
            Y = Y_tmp

            ! Save state every 'epochs' iterations
            if (mod(j,epochs) == 0) then
                call write_to_file(Y, t_min + dt*real((j), pr), dx, conversion_factors, unitnum)
            end if
        end do

        ! Save the last iteration if it has not been already saved
        if (write_last_iter) then
            call write_to_file(Y, t_min + dt*real((int((t_max - t_min)/dt)), pr), dx, conversion_factors, unitnum)
        end if
    close(unitnum)

end subroutine crank_nicolson_Neumann

subroutine analitical_sol(initial_conditions, Y, dt, t_min, t_max, filename, conversion_factors, epochs)
    real (kind=pr), dimension (:), intent (inout)               :: initial_conditions
    real (kind=pr), dimension (:), allocatable, intent (out)    :: Y
    real (kind=pr), dimension (:), intent(in)                   :: conversion_factors
    real (kind=pr), intent (in)                                 :: dt, t_max,t_min
    real (kind=pr), dimension (size(initial_conditions))        :: length_points
    integer (kind=int_medium), intent(in)                       :: epochs
    real (kind=pr)                                              :: dx
    integer (kind=int_huge)                                     :: j
    integer (kind=int_medium)                                   :: unitnum
    character (len=:), allocatable, intent(in)                  :: filename
    logical                                                     :: write_last_iter

    dx = 1._pr/real(size(initial_conditions)-1,pr)   ! dx = 1/(grid points -1)
    length_points = (/(dx*j, j = 0, size(initial_conditions)-1)/)

    ! Check if all iterations will be written or the last one needs to be saved separately
    if(mod(int((t_max - t_min)/dt),epochs)==0) then
        write_last_iter = .false.
    else
        write_last_iter = .true.
    end if

    open(newunit=unitnum, file=filename, status='replace')
        write(unitnum, *) "# time | x | Temperature "
        call write_to_file(initial_conditions, t_min, dx, conversion_factors, unitnum)
        do j = 1, int((t_max - t_min)/dt)
            Y  = analitical_solution(length_points , t_max)

            ! Save state every 'epochs' iterations
            if (mod(j,epochs) == 0) then
                call write_to_file(Y, t_min + dt*real((j), pr), dx, conversion_factors, unitnum)
            end if
        end do

        ! Save the last iteration if it has not been already saved
        if (write_last_iter) then
            call write_to_file(Y, t_min + dt*real((int((t_max - t_min)/dt)), pr), dx, conversion_factors, unitnum)
        end if
    close(unitnum)

end subroutine analitical_sol

END MODULE