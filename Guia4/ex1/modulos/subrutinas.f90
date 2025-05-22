MODULE subrutinas
use precision
use formats
use omp_lib
use mzranmod_threadsafe
implicit none
    integer(int_large)                  :: x_size, y_size
    type(MZRanState), allocatable       :: states(:)
    integer                             :: seeds(8,4)
    real(pr)                            :: N_spins
    integer                             :: autocorrelation_len_max
    integer                             :: i_last, i_next
    real(pr), allocatable               :: energy_buffer(:), magnetization_buffer(:)
    real(pr), allocatable               :: energy_autocorr(:)
    real(pr), allocatable               :: magnetization_autocorr(:)


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

subroutine lattice_init(lattice, initial_magnetization, state)
    integer, intent(out)     :: lattice(:,:)
    real(pr), intent(in)     :: initial_magnetization
    type(MZRanState)         :: state
    integer(int_large)       :: i, j

    do i = 1, x_size
        do j = 1, y_size
            if (rmzran_threadsafe(state) <= initial_magnetization) then
                lattice(i,j) = 1
            else
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

subroutine MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, transition_probability, state)
    integer, intent(inout), allocatable  :: lattice(:,:)
    real(pr), intent(inout)              :: Energy, magnetization
    real(pr), intent(in)                 :: transition_probability(-2:2)
    integer                              :: quarter_dE
    integer                              :: i, j, k, up, down, right, left
    type(MZRanState)                     :: state

    do k = 1, int(N_spins, int_large)

        ! Getting a random index within the bounds of the lattice indexes
        i = min(int(rmzran_threadsafe(state)*x_size + 1), x_size) ! Alternatively could use floor()
        j = min(int(rmzran_threadsafe(state)*y_size + 1), y_size)

        up    = mod(i, x_size) + 1
        down  = mod(i - 2 + x_size, x_size) + 1
        right = mod(j, y_size) + 1
        left  = mod(j - 2 + y_size, y_size) + 1

        ! The energy difference will only depend on the nearest neightbours' interaction with the flipped spin as follows:
        quarter_dE = lattice(i,j) * (lattice(up,j) + lattice(down,j) + lattice(i,right) + lattice(i,left))/2


        if (quarter_dE<=0) then
            lattice(i,j) = -lattice(i,j)
            magnetization = magnetization + real(2*lattice(i,j),pr)
            Energy = Energy + real(4*quarter_dE,pr)
        else if (rmzran_threadsafe(state) < transition_probability(quarter_dE)) then ! Note that this way avoid extra unnecessary calculations
            lattice(i,j) = -lattice(i,j)
            magnetization = magnetization + real(2*lattice(i,j),pr)
            Energy = Energy + real(4*quarter_dE,pr)
        end if
    end do

end subroutine MonteCarlo_step_PARALLEL

subroutine update_observables_absMagnetization(energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle &
    , magnetization_per_particle)
    real(pr)                :: Energy, magnetization
    real(kind=pr)           :: u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle, magnetization_per_particle

    energy_per_particle = Energy/N_spins
    magnetization_per_particle = abs(magnetization)/N_spins

    u_avg = u_avg + energy_per_particle
    uSqr_avg = uSqr_avg + energy_per_particle*energy_per_particle

    m_avg = m_avg + magnetization_per_particle
    mSqr_avg = mSqr_avg + magnetization_per_particle*magnetization_per_particle

end subroutine update_observables_absMagnetization

subroutine update_observables_normalMagnetization(energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle &
, magnetization_per_particle)
    real(pr)                :: Energy, magnetization
    real(kind=pr)           :: u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle, magnetization_per_particle

    magnetization_per_particle = magnetization/N_spins
    energy_per_particle = Energy/N_spins
    u_avg = u_avg + energy_per_particle
    uSqr_avg = uSqr_avg + energy_per_particle*energy_per_particle
    m_avg = m_avg + magnetization_per_particle
    mSqr_avg = mSqr_avg + magnetization_per_particle*magnetization_per_particle

end subroutine update_observables_normalMagnetization

subroutine init_autocorr()
    allocate(energy_buffer(0:autocorrelation_len_max))
    allocate(magnetization_buffer(0:autocorrelation_len_max))
    allocate(energy_autocorr(0:autocorrelation_len_max))
    allocate(magnetization_autocorr(0:autocorrelation_len_max))
    energy_buffer = 0.0_pr
    magnetization_buffer = 0.0_pr
    energy_autocorr = 0.0_pr
    magnetization_autocorr = 0.0_pr
end subroutine init_autocorr

subroutine update_autocorrelation_contributions(magnetization_per_particle, energy_per_particle, energy_autocorr &
    , magnetization_autocorr, autocorr_count)
    real(pr), intent(in)       :: magnetization_per_particle, energy_per_particle
    real(pr), intent(inout)    :: energy_autocorr(0:autocorrelation_len_max)
    real(pr), intent(inout)    :: magnetization_autocorr(0:autocorrelation_len_max)
    integer                    :: j
    integer, intent(inout)     :: autocorr_count

    autocorr_count = autocorr_count + 1
    i_last = mod(autocorr_count, autocorrelation_len_max) + 1
    energy_buffer(i_last) = energy_per_particle
    magnetization_buffer(i_last) = magnetization_per_particle

    !if (autocorr_count >= autocorrelation_len_max) then
        do j = 0, min(autocorrelation_len_max , autocorr_count-1)
            i_next = mod(autocorr_count - j, autocorrelation_len_max) + 1
            energy_autocorr(j) = energy_autocorr(j) + energy_buffer(i_last) * energy_buffer(i_next)
            magnetization_autocorr(j) = magnetization_autocorr(j) + magnetization_buffer(i_last)*magnetization_buffer(i_next)
        end do
   ! end if

end subroutine update_autocorrelation_contributions

subroutine save_autocorrelation(energy_autocorr, magnetization_autocorr, KbT, threadID)
    real(pr), intent(in)    :: energy_autocorr(0:autocorrelation_len_max)
    real(pr), intent(in)    :: magnetization_autocorr(0:autocorrelation_len_max)
    integer(int_small)      :: threadID
    integer                 :: unit_autocorrelation, i
    character(len=34)       :: file_autocorrelation
    real(pr)                :: KbT

    call create_file_name("datos/autocorrelation_T_", KbT, ".out", file_autocorrelation)
    open(newunit=unit_autocorrelation, file=file_autocorrelation, status='replace')
        write(unit_autocorrelation,format_style_header) "##  Cell dimensions: ", x_size,"x",y_size
        write(unit_autocorrelation,format_style_header) "## Thread ID = ", threadID
        write(unit_autocorrelation,'(a)') "## autocorrelation length | magnetization autocorrelation | energy "// &
            "autocorrelation"
        do i = 0, autocorrelation_len_max
            write(unit_autocorrelation,format_style1) i, magnetization_autocorr(i), energy_autocorr(i)
        end do
    close(unit_autocorrelation)

end subroutine save_autocorrelation

subroutine append_lattice_binary(filename, lattice)
    character(len=*), intent(in)    :: filename
    integer, intent(in)             :: lattice(:,:)
    integer                         :: i, j, k, bit_index
    integer                         :: unit
    integer(kind=1), allocatable    :: byte_array(:)
    integer                         :: num_bits, num_bytes

    num_bits = x_size * y_size
    num_bytes = (num_bits + 7) / 8  ! Ensure to get the right number of bytes
    allocate(byte_array(num_bytes))
    byte_array = 0

    k = 1
    bit_index = 0
    do j = 1, y_size
        do i = 1, x_size
            if (lattice(i, j) == 1) then
                byte_array(k) = ibset(byte_array(k), mod(bit_index, 8))   ! flip the corresponding bit to 1
            end if
            bit_index = bit_index + 1
            if (mod(bit_index, 8) == 0) k = k + 1
        end do
    end do

    open(newunit=unit, file=filename, form='unformatted', access='stream', &
         status='old', action='write', position='append')
    write(unit) byte_array
    close(unit)
end subroutine

END MODULE