!***************************************************************
! Program: ex1.f90
! Purpose: Analyzing the Ising model and different observables as a function of the dimensionless temperature Kb*T/J.
!
!
! Description:  The following abbreviations for different variables' names is implements: Sqr = square, avg = average, var = variance, inv = inverse
!
! Input: The input information is specified in an "input.nml" file. Default values are provided for every variable in the "Input settings" section, and an example input file is also provided.
!
!
! Output: all output data is written to a "datos" folder which MUST exist before the calculations. Every file has the suffix ".out". Observables and errors are printed in "datos/temperature_functions.out" as they are calculated, and ordered in "datos/temperature_functions_sorted.out" after the calculations have finished. Furthermore, the associated errors are also stored, but not the standard deviations, as a simple multiplication by sqrt(N-1) suffices to calculate it if necessary in post processing. The thermalization data is saved (if specified in the input.nml) to files specifying the temperature used and the initial magnetization chosen.
!
!
!
! Room for improvement: The correlation of the seeds chosen must be studied for the reliability of the calculations. The seeds used are simply an example and by no means the best choice.
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precision
    use subrutinas
    use observables
    use omp_lib
    use parsing
    use mzranmod_threadsafe
    implicit none

    real(kind=pr)                   :: u_avg, uSqr_avg, u_var, u_error, m_avg, mSqr_avg, m_var, m_error, Binder
    real(kind=pr)                   :: energy_per_particle, magnetization_per_particle, real_MC_steps
    real(kind=pr)                   :: susceptibility, capacity
    real(kind=pr), allocatable      :: beta(:), KbT(:)
    real(pr)                        :: transition_probability(-2:2)
    integer, allocatable            :: lattice(:,:)
    integer(int_small)              :: threadID
    integer                         :: i, j, k, l, unit_steps, unit_temperature, nthreads, status, frame
    real(pr)                        :: Energy, magnetization
    integer(int_large)              :: num_measurements, autocorr_count
    character(len=31)               :: file_temperature
    character(len=28)               :: file_steps, file_lattice
    character(len=8)                :: prefix
    character(len=14)               :: suffix
    character(len=140)              :: sort_command


    abstract interface
        subroutine update(E, M, u_avg, u_var, m_avg, m_var, EpP, MpP, Binder)
            use precision
            real(pr), intent(in)    :: E, M
            real(pr), intent(inout) :: u_avg, u_var, m_avg, m_var, EpP, MpP, Binder
        end subroutine update
    end interface
    procedure (update), pointer :: update_observables => null()

!###################################################################################################
!   Set default values and parse input file
!###################################################################################################

    call set_defaults()
    call parse_input()

!###################################################################################################
!   Initializing other relevant parameters
!###################################################################################################

    ! Random number generator
    call MZRanState_init()
    nthreads = size(states)

    ! Cell parameters
    allocate(lattice(x_size,y_size))
    N_spins = real(x_size*y_size,pr)

    ! Temperatures to examine
    select case (T_range)
        case(.true.)
            allocate(KbT(KbT_steps))
            KbT = (/(Kbt_min + real(i,pr)*(KbT_max - KbT_min)/real(KbT_steps-1,pr), i=0, KbT_steps-1)/)
        case(.false.)
            allocate(KbT(1))
            KbT = KbT_user
    end select
    beta = 1._pr/KbT

    ! Net MC steps
    num_measurements = (MC_steps - transitory_steps) / step_jump
    real_MC_steps = real(num_measurements, pr)

    ! Select whether the absolute value of the magnetization is to be averaged
    select case (use_absolute_magnetization)
        case(.true.)
            if (do_binder) update_observables => update_observables_andBinder_absMagnetization
            if (.not. do_binder) update_observables => update_observables_absMagnetization
        case(.false.)
            if (do_binder) update_observables => update_observables_andBinder_normalMagnetization
            if (.not. do_binder) update_observables => update_observables_normalMagnetization
    end select

    ! Only save lattice frames if one temperature is being inspected
    if (T_range) save_lattice_evolution = .false.

    ! Files' names and parts
    prefix = "datos/T_"
    call create_suffix("_m0_", initial_magnetization, ".out", suffix)
    file_temperature = "datos/temperature_functions.out"
    if (do_autocorrelation)  call init_autocorr()

!##################################################################################
!    Begin and save calculations for all T specified in the KbT array
!##################################################################################

    open(newunit=unit_temperature, file=file_temperature, status='unknown')
    write(unit_temperature,format_style_header) "##  Cell dimensions: ", x_size,"x",y_size
    write(unit_temperature,'(a,I7)') "##  Effective Monte Carlo steps:", num_measurements
    if (do_binder) then
        write(unit_temperature,'(a)') "##     KbT      | Thread ID |       <m>      |     m Error    | susceptibility |"//&
        "      <u>      |    u Error     |    capacity     |  Binder cumulant"
    else
          write(unit_temperature,'(a)') "##     KbT      | Thread ID |       <m>      |     m Error    | susceptibility |"//&
        "      <u>      |    u Error     |    capacity"
    end if

    !$omp parallel do private(Energy, magnetization, lattice, u_avg, uSqr_avg, m_avg, mSqr_avg, threadID, transition_probability &
    !$omp   , j, i, k, l, magnetization_per_particle, energy_per_particle, unit_steps, file_steps, energy_buffer, Binder &
    !$omp   , magnetization_buffer, autocorr_count, i_last, i_next, energy_autocorr, magnetization_autocorr, file_lattice, frame) &
    !$omp shared(KbT, nthreads, unit_temperature, format_style0, format_style1, format_style_header, beta &
    !$omp   , states, seeds, real_MC_steps, prefix, suffix, N_spins, step_jump, autocorrelation_len_max &
    !$omp   , save_lattice_evolution, lattice_frames) &
    !$omp   schedule(dynamic)


    do j = 1, size(KbT)

        ! Initialize RNG states with unique seeds
        threadID = int(omp_get_thread_num() + 1, int_small)
        call init_mzran_threadsafe(states(threadID), seeds(threadID,1), seeds(threadID,2), seeds(threadID,3), seeds(threadID,4))

        ! Initialize system and variables
        call lattice_init(lattice, initial_magnetization, states(threadID))
        call get_lattice_energy_vectorized(lattice, Energy)
        magnetization = sum(real(lattice,pr))
        u_avg = 0._pr
        uSqr_avg = 0._pr
        m_avg = 0._pr
        mSqr_avg = 0._pr
        transition_probability = (/(exp (-beta(j)*real(4*i,pr)), i=-2,2)/)

        ! Select wether to save the steps to files or not
        select case (save_thermalization)
            case(.true.)
                call create_file_name(prefix, KbT(j), suffix, file_steps)
                print*, "Calculating KbT = ",  KbT(j), "Filename:  ", file_steps, " ThreadID:", threadID
                open(newunit=unit_steps, file=file_steps, status='replace')
                    write(unit_temperature,format_style_header) "##  Cell dimensions: ", x_size,"x",y_size

                    magnetization_per_particle = magnetization/N_spins
                    energy_per_particle = energy/N_spins
                    write(unit_steps,format_style_header) "## Thread ID = ", threadID
                    write(unit_steps,'(a)') "## MC steps | energy per particle | magnetization per particle"
                    write(unit_steps,format_style1) 0, energy_per_particle, magnetization_per_particle

                    ! Transitory steps
                    do i = 1, transitory_steps
                        call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, transition_probability, states(threadID))
                        magnetization_per_particle = magnetization/N_spins
                        energy_per_particle = Energy/N_spins
                        if (save_transitory) write(unit_steps,format_style1) i, energy_per_particle, magnetization_per_particle
                    end do

                    ! Initialize autocorrelation variables
                    if (do_autocorrelation)  then
                        energy_buffer = 0.0_pr
                        magnetization_buffer = 0.0_pr
                        energy_autocorr = 0.0_pr
                        magnetization_autocorr = 0.0_pr
                        autocorr_count = 0
                    end if

                    ! Create file to save lattice states and initialize frame counter
                    if (save_lattice_evolution) then
                        call create_file_name("datos/lattice_T_", KbT(j), ".bin", file_lattice)
                        open(newunit=i, file=file_lattice, status="replace", action="write")
                        close(i)
                        frame = 0
                    end if

                    ! Important steps
                    do i = 1, num_measurements
                        do k = 1, step_jump
                            call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, transition_probability, states(threadID))
                        end do

                        call update_observables(Energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg &
                        , energy_per_particle, magnetization_per_particle, Binder)

                        write(unit_steps,format_style1) i, energy_per_particle, magnetization_per_particle
                        if (do_autocorrelation)  call update_autocorrelation_contributions(magnetization_per_particle &
                            , energy_per_particle, energy_autocorr, magnetization_autocorr, autocorr_count)

                        if (save_lattice_evolution) then
                            call append_lattice_binary(file_lattice, lattice)
                            frame = frame + 1
                            if (frame == lattice_frames) save_lattice_evolution = .false.
                        end if
                end do
                close(unit_steps)

            case(.false.)
                print*, "Calculating KbT = ",  KbT(j), " ThreadID:", threadID

                ! Transitory steps
                do i = 1, transitory_steps
                    call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, transition_probability, states(threadID))
                end do

                ! Initialize autocorrelation variables
                if (do_autocorrelation)  then
                    energy_buffer = 0.0_pr
                    magnetization_buffer = 0.0_pr
                    energy_autocorr = 0.0_pr
                    magnetization_autocorr = 0.0_pr
                    autocorr_count = 0
                end if

                ! Create file to save lattice states and initialize frame counter
                if (save_lattice_evolution) then
                    call create_file_name("datos/lattice_T_", KbT(j), ".bin", file_lattice)
                    open(newunit=i, file=file_lattice, status="replace", action="write")
                    close(i)
                    frame = 0
                end if

                ! Important steps
                do i = 1, num_measurements
                    do k = 1, step_jump
                        call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, transition_probability, states(threadID))
                    end do

                    call update_observables(Energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle&
                        , magnetization_per_particle, Binder)

                    if (do_autocorrelation)  call update_autocorrelation_contributions(magnetization_per_particle &
                        , energy_per_particle, energy_autocorr, magnetization_autocorr, autocorr_count)

                    if (save_lattice_evolution) then
                        call append_lattice_binary(file_lattice, lattice)
                        frame = frame + 1
                        if (frame == lattice_frames) save_lattice_evolution = .false.
                    end if
                end do
        end select

        ! Calculate statistics for u and m
        u_avg = u_avg / real_MC_steps
        uSqr_avg = uSqr_avg / real_MC_steps
        u_var = (uSqr_avg - u_avg*u_avg)
        u_error = sqrt(u_var/(real_MC_steps-1_pr))

        m_avg = m_avg / real_MC_steps
        mSqr_avg = mSqr_avg / real_MC_steps
        m_var = (mSqr_avg - m_avg*m_avg)
        m_error = sqrt(m_var/(real_MC_steps-1_pr))
        if (do_binder) Binder = 1._pr - Binder/(real_MC_steps*3._pr*mSqr_avg*mSqr_avg)

        ! Calculate other observables
        capacity = N_spins*u_var*beta(j)*beta(j)
        susceptibility = N_spins*m_var*beta(j)

        !$omp critical
            if (do_binder) then
                write(unit_temperature,format_style2) KbT(j), threadID, m_avg, m_error, susceptibility, u_avg, u_error, capacity &
                    , Binder
            else
                write(unit_temperature,format_style2) KbT(j), threadID, m_avg, m_error, susceptibility, u_avg, u_error, capacity
            end if
            call flush(unit_temperature)
        !$omp end critical

        if (do_autocorrelation) then
            ! Subtract mean and divide by variance and normalize to get autocorrelation A(k) for energy and magnetization per particle
            energy_autocorr = (/(energy_autocorr(i) / real(autocorr_count + 1 - i, pr), i=0 &
                , autocorrelation_len_max)/)
            magnetization_autocorr = (/(magnetization_autocorr(i) / real(autocorr_count + 1 - i, pr), i=0 &
                , autocorrelation_len_max)/)
            energy_autocorr = (energy_autocorr - u_avg*u_avg) / u_var
            magnetization_autocorr = (magnetization_autocorr - m_avg*m_avg) / m_var

            call save_autocorrelation(energy_autocorr, magnetization_autocorr, KbT(j), threadID)

        end if

    end do

    !$omp end parallel do

    sort_command = "(head -n 1 datos/temperature_functions.out && tail -n +2 datos/temperature_functions.out | sort -g) "//&
   "> datos/temperature_functions_sorted.out"
    call execute_command_line(sort_command, wait=.true., exitstat=status)
    close(unit_temperature)
    deallocate(KbT)


end program ex1