MODULE subrutinas
use precision
use formats
use omp_lib
use mzranmod
use mzranmod_threadsafe
implicit none
    integer(int_large)                  :: x_size, y_size
    type(MZRanState), allocatable       :: states(:)
    integer                             :: seeds(8,4)
    real(pr)                            :: N_spinors


contains

subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                :: fmt  ! Format descriptor
    character(len=12)               :: x1   ! Temporary string for formatted real
    character(len=*), intent(in)    :: prefix, suffix
    character(len=*)                :: filename
    real(kind=pr), intent(in)       :: num  ! Input real number

    fmt = '(F8.4)'  ! Adjust as needed
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

subroutine create_suffix(prefix, num, suffix, suffix_out)
    character(len=8)                :: fmt  ! Format descriptor
    character(len=12)               :: x1   ! Temporary string for formatted real
    character(len=*), intent(in)    :: prefix, suffix
    character(len=*)                :: suffix_out
    real(kind=pr), intent(in)       :: num  ! Input real number

    fmt = '(F6.4)'  ! Adjust as needed
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    suffix_out = prefix // trim(x1) // suffix

end subroutine create_suffix

subroutine MZRanState_init()
    integer      :: nthreads
    ! Determine number of threads and allocate state array
    !$omp parallel
    !$omp single
    nthreads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    allocate(states(nthreads))

    seeds(1,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(2,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(3,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(4,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(5,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(6,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(7,:) = (/521288629, 362436069, 16163801, 1131199299/)
    seeds(8,:) = (/521288629, 362436069, 16163801, 1131199299/)

end subroutine MZRanState_init

subroutine lattice_init(lattice, initial_magnetization)
    integer, intent(out)     :: lattice(:,:)
    real(pr), intent(in)                :: initial_magnetization
    real(kind=pr)                       :: rnd_num
    integer(int_large)                  :: i, j

    do i = 1, x_size
        do j = 1, y_size
            rnd_num = rmzran()
            if (rnd_num <= initial_magnetization) then
                lattice(i,j) = 1
            elseif (rnd_num > initial_magnetization) then
                lattice(i,j) = -1
            end if
        end do
    end do

end subroutine lattice_init

subroutine get_lattice_energy(lattice, Energy)
    integer, intent(in), allocatable     :: lattice(:,:)
    real(pr)                             :: Energy
    integer(int_large)                   :: i, j, i_up, j_right !, j_left, i_down

    Energy = 0

    do i = 1, x_size
        do j = 1, y_size
            ! Periodic boundary conditions
            i_up    = mod(i, x_size) + 1
            !i_down  = mod(i - 2 + x_size, x_size) + 1
            j_right = mod(j, y_size) + 1
            !j_left  = mod(j - 2 + y_size, y_size) + 1

            Energy = Energy - real(lattice(i,j) * (lattice(i_up,j) + lattice(i,j_right)),pr) !+ lattice(i,j_left)+ lattice(i_down,j) )
        end do
    end do

    !Energy = Energy/2  ! If we overcount, E is even and thus the division has no remainder

end subroutine get_lattice_energy

subroutine get_lattice_energy_vectorized(lattice, Energy)
    integer, intent(in), allocatable                :: lattice(:,:)
    real(pr)                                        :: Energy
    integer(int_large)                              :: i
    integer(int_large), dimension(x_size,y_size)    :: shifted_right, shifted_left, shifted_up, shifted_down

    ! Shift the lattice using periodic boundary conditions to vectorize calculations
    shifted_right = lattice(mod( (/(i, i=1,x_size)/) ,x_size ) + 1, :)      ! Shift right
    shifted_left = lattice(mod((/(i-2+x_size, i=1,x_size)/),x_size) + 1, :)   ! Shift left
    shifted_up = lattice(:, mod((/(i, i=1,y_size)/), y_size) + 1)      ! Shift up
    shifted_down = lattice(:, mod((/(i-2+y_size, i=1,y_size)/), y_size) + 1)  ! Shift down

    Energy = -SUM(int(lattice * (shifted_right + shifted_left + shifted_up + shifted_down),int_large)) /2 ! If we overcount, E is even and thus the division has no remainder

end subroutine get_lattice_energy_vectorized

subroutine MonteCarlo_step(lattice, Energy, magnetization, beta)
    integer, intent(inout), allocatable  :: lattice(:,:)
    real(pr), intent(inout)              :: Energy, magnetization
    real(pr), intent(in)                 :: beta
    integer                              :: dE
    real(kind=pr)                        :: threshold
    integer                              :: i, j, k, up, down, right, left

    do k = 1, int(N_spinors, int_large)

        ! Getting a random index within the bounds of the lattice indexes
        i = min(int(rmzran()*x_size + 1), x_size) ! Alternatively could use floor()
        j = min(int(rmzran()*y_size + 1), y_size)

        up    = mod(i, x_size) + 1
        down  = mod(i - 2 + x_size, x_size) + 1
        right = mod(j, y_size) + 1
        left  = mod(j - 2 + y_size, y_size) + 1

        ! The energy difference will only depend on the nearest neightbours' interaction with the flipped spin as follows:
        dE = 2*lattice(i,j) * (lattice(up,j) + lattice(down,j) + lattice(i,right) + lattice(i,left))

        threshold = exp (-beta*real(dE,pr))
        if (dE<=0) then
            lattice(i,j) = -lattice(i,j)
            magnetization = magnetization + real(2*lattice(i,j),pr)
            Energy = Energy + real(dE,pr)
        else if (rmzran() < threshold) then
            lattice(i,j) = -lattice(i,j)
            magnetization = magnetization + real(2*lattice(i,j),pr)
            Energy = Energy + real(dE,pr)
        end if
    end do

end subroutine MonteCarlo_step

subroutine MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta, state)
    integer, intent(inout), allocatable  :: lattice(:,:)
    real(pr), intent(inout)              :: Energy, magnetization
    real(pr), intent(in)                 :: beta
    integer                              :: dE
    real(kind=pr)                        :: threshold
    integer                              :: i, j, k, up, down, right, left
    type(MZRanState)                     :: state

    do k = 1, int(N_spinors, int_large)

        ! Getting a random index within the bounds of the lattice indexes
        i = min(int(rmzran_threadsafe(state)*x_size + 1), x_size) ! Alternatively could use floor()
        j = min(int(rmzran_threadsafe(state)*y_size + 1), y_size)

        up    = mod(i, x_size) + 1
        down  = mod(i - 2 + x_size, x_size) + 1
        right = mod(j, y_size) + 1
        left  = mod(j - 2 + y_size, y_size) + 1

        ! The energy difference will only depend on the nearest neightbours' interaction with the flipped spin as follows:
        dE = 2*lattice(i,j) * (lattice(up,j) + lattice(down,j) + lattice(i,right) + lattice(i,left))

        threshold = exp (-beta*real(dE,pr))
        if (dE<=0) then
            lattice(i,j) = -lattice(i,j)
            magnetization = magnetization + real(2*lattice(i,j),pr)
            Energy = Energy + real(dE,pr)
        else if (rmzran_threadsafe(state) < threshold) then
            lattice(i,j) = -lattice(i,j)
            magnetization = magnetization + real(2*lattice(i,j),pr)
            Energy = Energy + real(dE,pr)
        end if
    end do

end subroutine MonteCarlo_step_PARALLEL

subroutine update_observables_absMagnetization(energy, magnetization, u_avg, u_var, m_avg, m_var, energy_per_particle &
, magnetization_per_particle)
    real(pr)                :: Energy, magnetization
    real(kind=pr)           :: u_avg, u_var, m_avg, m_var, energy_per_particle, magnetization_per_particle

    magnetization_per_particle = magnetization/N_spinors
    energy_per_particle = Energy/N_spinors
    u_avg = u_avg + energy_per_particle
    u_var = u_var + energy_per_particle*energy_per_particle
    m_avg = m_avg + abs(magnetization_per_particle)
    m_var = m_var + magnetization_per_particle*magnetization_per_particle

end subroutine update_observables_absMagnetization

subroutine update_observables_normalMagnetization(energy, magnetization, u_avg, u_var, m_avg, m_var, energy_per_particle &
, magnetization_per_particle)
    real(pr)                :: Energy, magnetization
    real(kind=pr)           :: u_avg, u_var, m_avg, m_var, energy_per_particle, magnetization_per_particle

    magnetization_per_particle = magnetization/N_spinors
    energy_per_particle = Energy/N_spinors
    u_avg = u_avg + energy_per_particle
    u_var = u_var + energy_per_particle*energy_per_particle
    m_avg = m_avg + abs(magnetization_per_particle)
    m_var = m_var + magnetization_per_particle*magnetization_per_particle

end subroutine update_observables_normalMagnetization

END MODULE