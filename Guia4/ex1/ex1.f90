!***************************************************************
! Program: ex1.f90
! Purpose: Analysing the Ising model and different observables as a function of the dimentionaless temperature Kb*T/J.
!
!
! Description:  The following abreviations for different variables' names is implementes: Sqr = square, avg = average, var = variance, inv = inverse
!
! Input: The imput information is specified in an "input.nml" file. Defauls values are provided for every variable in the "Input settings" section, and an example input file is also provided.
!
!
! Output: all output data is written to a "datos" folder which MUST exist berofe the calculations. Every file has the suffix ".out". Observables and errors are printed in "datos/temperature_functions.out" as they are calculated, and ordered in "datos/temperature_functions_sorted.out" after the calculations have finished. The thermalization data is saved (if specified in the input.nml) to files specifying the temperature used and the initial magnetization chosen.
!
!
!
! Room for improvement: The correlation of the seeds chones must be studied for the reliability of the calculations. The seeds used are simply an example and by no means the best choise.
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precision
    use subrutinas
    use omp_lib
    use mzranmod_threadsafe
    use mzranmod
    implicit none
    real(kind=pr)                   :: u_avg, uSqr_avg, u_var, u_error, m_avg, mSqr_avg, m_var, m_error
    real(kind=pr)                   :: energy_per_particle, magnetization_per_particle, real_MC_steps
    real(kind=pr)                   :: susceptibility, capacity, KbT_min, KbT_max, KbT_user, initial_magetization, N_spinors
    real(kind=pr), allocatable      :: beta(:), KbT(:)
    integer, allocatable            :: lattice(:,:)
    integer(int_small)              :: threadID
    integer                         :: i, j, unit_steps, unit_temperature, unitnum, nthreads, status
    real(pr)                        :: Energy, magnetization
    integer(int_large)              :: MC_steps, transitory_steps, KbT_steps
    character(len=31)               :: file_temperature
    character(len=28)               :: file_steps
    character(len=8)                :: prefix
    character(len=14)               :: suffix
    character(len=140)              :: command
    logical                         :: T_range, save_thermalization, use_abolute_magnetization

    abstract interface
        subroutine update(E, M, N, u_avg, u_var, m_avg, m_var, EpP, MpP)
            use precision
            real(pr)                :: E, M
            real(pr)                :: N
            real(kind=pr)           :: u_avg, u_var, m_avg, m_var, EpP, MpP
        end subroutine update
    end interface
    procedure (update), pointer :: update_observables => null()

!##################################################################################################
!       Input settings
!##################################################################################################

    ! Namelist blocks
    namelist /physical/ KbT_min, KbT_max, KbT_steps, T_range, KbT_user, initial_magetization, x_size, y_size
    namelist /calculation/ MC_steps, transitory_steps, save_thermalization, use_abolute_magnetization

    ! DEFAULT SETTINGS
        ! Physical problems' characteristics
            KbT_user    = 2.2676_pr
            T_range     = .false.
            KbT_min     = 0.1_pr
            KbT_max     = 3.5
            KbT_steps   = 33
            x_size      = 40
            y_size      = 40
            initial_magetization = 1._pr

        !Calculation settings
            MC_steps = 1000000
            transitory_steps = MC_steps/2
            save_thermalization = .false.
            use_abolute_magnetization = .true.

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=physical)
        read(unitnum, nml=calculation)
    close(unitnum)

!###################################################################################################
!   Initializing other relevant parameters
!###################################################################################################

    call MZRanState_init() ! Initialize the "states" array for the parallelization
    nthreads = size(states)

    allocate(lattice(x_size,y_size))
    N_spinors = real(x_size*y_size,pr)

    select case (T_range)
        case(.true.)
            allocate(KbT(KbT_steps))
            KbT = (/(Kbt_min + real(i,pr)*(KbT_max - KbT_min)/real(KbT_steps-1,pr), i=0, KbT_steps-1)/)
        case(.false.)
            allocate(KbT(1))
            KbT = KbT_user
    end select
    beta = 1._pr/KbT
    real_MC_steps = real(MC_steps - transitory_steps,pr)

    select case (use_abolute_magnetization)
        case(.true.)
            update_observables => update_observables_absMagnetization
        case(.false.)
            update_observables => update_observables_normalMagnetization
    end select

    file_temperature = "datos/temperature_functions.out"


!##################################################################################
!    Begin and save calculations for all T specified in the KbT array
!##################################################################################

    open(newunit=unit_temperature, file=file_temperature, status='unknown')
    write(unit_temperature,*) "##     KbT      | Thread ID |       <m>      |     m Error    | susceptibility |"//&
    "      <u>      |    u Error     |    capacity"

    !$omp parallel do private(Energy, magnetization, lattice, u_avg, uSqr_avg, m_avg, mSqr_avg, &
    !$omp   magnetization_per_particle, energy_per_particle, unit_steps, file_steps, prefix, &
    !$omp   suffix, j, i, threadID) shared(KbT, nthreads, unit_temperature, format_style0, &
    !$omp   format_style1, beta, states, seeds, real_MC_steps) schedule(dynamic)

    do j = 1, size(KbT)
        ! This initialization should be done outside the loop but got "-Wuninitialized" variables
        prefix = "datos/T_"
        call create_suffix("_m0_", initial_magetization, ".out", suffix)
    !    suffix = "_T0_zero.out"


        ! Initialize RNG states with unique seeds
        threadID = int(omp_get_thread_num() + 1, int_small)
        call init_mzran_threadsafe(states(threadID), seeds(threadID,1), seeds(threadID,2), seeds(threadID,3), seeds(threadID,4))

        call lattice_init(lattice, initial_magetization)
        call get_lattice_energy_vectorized(lattice, Energy)
        magnetization = sum(int(lattice,int_large))
        u_avg = 0._pr
        uSqr_avg = 0._pr
        m_avg = 0._pr
        mSqr_avg = 0._pr



        select case (save_thermalization)
            case(.true.)
                call create_file_name(prefix, KbT(j), suffix, file_steps)
                print*, "Calculating KbT = ",  KbT(j), "Filename:  ", file_steps, " ThreadID:", threadID
                open(newunit=unit_steps, file=file_steps, status='replace')

                    magnetization_per_particle = real(magnetization,pr)/N_spinors
                    energy_per_particle = real(energy,pr)/N_spinors
                    write(unit_steps,*) "## Thread ID = ", threadID
                    write(unit_steps,*) "## MC steps | energy per particle | magnetization per particle"
                    write(unit_steps,format_style1) 0, energy_per_particle, magnetization_per_particle

                    do i = 1, transitory_steps
                        call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta(j), states(threadID))
                        magnetization_per_particle = real(magnetization,pr)/N_spinors
                        energy_per_particle = real(energy,pr)/N_spinors
                        write(unit_steps,format_style1) i, energy_per_particle, magnetization_per_particle
                    end do

                    do i = transitory_steps + 1, MC_steps
                        call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta(j), states(threadID))
                        call update_observables(Energy, magnetization, N_spinors,u_avg, uSqr_avg, m_avg, mSqr_avg &
                        , energy_per_particle, magnetization_per_particle)
                        write(unit_steps,format_style1) i, energy_per_particle, magnetization_per_particle
                    end do
                close(unit_steps)

            case(.false.)
                print*, "Calculating KbT = ",  KbT(j), " ThreadID:", threadID

                do i = 1, transitory_steps
                    call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta(j), states(threadID))
                end do

                do i = transitory_steps + 1, MC_steps
                    call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta(j), states(threadID))
                    call update_observables(Energy, magnetization, N_spinors,u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle&
                    , magnetization_per_particle)
                end do
        end select

        u_avg = u_avg / real_MC_steps
        uSqr_avg = uSqr_avg / real_MC_steps
        u_var = uSqr_avg - u_avg*u_avg

        m_avg = m_avg / real_MC_steps
        mSqr_avg = mSqr_avg / real_MC_steps
        m_var = mSqr_avg - m_avg*m_avg

        capacity = u_var*beta(j)*beta(j)
        susceptibility = m_var*beta(j)

        u_error = sqrt(u_var/(real_MC_steps-1_pr))
        m_error = sqrt(m_var/(real_MC_steps-1_pr))

        !$omp critical
            write(unit_temperature,format_style2) KbT(j), threadID, m_avg, m_error, susceptibility, u_avg, capacity, u_error
            call flush(unit_temperature)
        !$omp end critical

    end do

    !$omp end parallel do

    command = "(head -n 1 datos/temperature_functions.out && tail -n +2 datos/temperature_functions.out | sort -g) "//&
   "> datos/temperature_functions_sorted.out"
    call execute_command_line(command, wait=.true., exitstat=status)
    close(unit_temperature)
    deallocate(KbT)


end program ex1