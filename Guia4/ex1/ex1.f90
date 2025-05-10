!***************************************************************
! Program: ex1.f90
! Purpose:
!
!
! Description:
!
! Input:
!
!
! Output:
!
!
!
! Room for improvement: We can modify the program to calculate in parallel different thermal baths (KbT) and initial conditions
!   by increasing in 1 the dimension of the arrays and scalars used.
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
    real(kind=pr)                   :: u_avg, u_var, m_avg, m_var, energy_per_particle, magnetization_per_particle, real_MC_steps
    real(kind=pr)                   :: susceptibility, capacity, KbT_min, KbT_max, KbT_user, initial_magetization, N_spinors
    real(kind=pr), allocatable      :: beta(:), KbT(:)
    integer, allocatable            :: lattice(:,:)
    integer(int_small)              :: threadID
    integer                         :: i, j, unit_steps, unit_temperature, unitnum, nthreads, status
    integer(int_large)              :: Energy, magnetization, MC_steps, transitory_steps, KbT_steps
    character(len=31)               :: file_temperature
    character(len=28)               :: file_steps
    character(len=8)                :: prefix
    character(len=14)               :: suffix
    character(len=140)              :: command
    logical                         :: T_range, save_thermalization

!##################################################################################################
!       Input settings
!##################################################################################################

    ! Namelist blocks
    namelist /physical/ KbT_min, KbT_max, KbT_steps, T_range, KbT_user, initial_magetization, x_size, y_size
    namelist /calculation/ MC_steps, transitory_steps, save_thermalization

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
    if (T_range) then
        allocate(KbT(KbT_steps))
        KbT = (/(Kbt_min + real(i,pr)*(KbT_max - KbT_min)/real(KbT_steps-1,pr), i=0, KbT_steps-1)/)
    else
        allocate(KbT(1))
        KbT = KbT_user
    end if
    beta = 1._pr/KbT
    real_MC_steps = real(MC_steps - transitory_steps,pr)

    file_temperature = "datos/temperature_functions.out"


!##################################################################################
!    Begin and save calculations for all T specified in the KbT array
!##################################################################################

    open(newunit=unit_temperature, file=file_temperature, status='unknown')
    write(unit_temperature,*) "##     KbT      | Thread ID |       <m>          | susceptibility |      <u>     | capacity"

    !$omp parallel do private(Energy, magnetization, lattice, u_avg, u_var, m_avg, m_var, &
    !$omp   magnetization_per_particle, energy_per_particle, unit_steps, file_steps, prefix, &
    !$omp   suffix, j, i, threadID) shared(KbT, nthreads, unit_temperature, format_style0, &
    !$omp   format_style1, beta, states, seeds, real_MC_steps) schedule(dynamic)

    do j = 1, size(KbT)
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
        u_var = 0._pr
        m_avg = 0._pr
        m_var = 0._pr



        select case (save_thermalization)
            case(.true.)
                call create_file_name(prefix, KbT(j), suffix, file_steps)
                print*, "Calculating KbT = ",  KbT(j), "Filename:  ", file_steps, " ThreadID:", threadID
                open(newunit=unit_steps, file=file_steps, status='replace')

                    magnetization_per_particle = real(magnetization,pr)/N_spinors
                    energy_per_particle = real(energy,pr)/N_spinors
                    write(unit_steps,*) "## Thread ID = ", threadID," | MC steps | energy per particle | magnetization_per_particle"
                    write(unit_steps,format_style1) 0, energy_per_particle, magnetization_per_particle

                    do i = 1, transitory_steps
                        call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta(j), states(threadID))
                        magnetization_per_particle = real(magnetization,pr)/N_spinors
                        energy_per_particle = real(energy,pr)/N_spinors
                        write(unit_steps,format_style1) i, energy_per_particle, magnetization_per_particle
                    end do

                    do i = transitory_steps + 1, MC_steps
                        call MonteCarlo_step_PARALLEL(lattice, Energy, magnetization, beta(j), states(threadID))
                        call update_observables(Energy, magnetization, N_spinors,u_avg, u_var, m_avg, m_var, energy_per_particle&
                        , magnetization_per_particle)
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
                    call update_observables(Energy, magnetization, N_spinors,u_avg, u_var, m_avg, m_var, energy_per_particle&
                    , magnetization_per_particle)
                end do
        end select

        u_avg = u_avg / real_MC_steps
        u_var = u_var / real_MC_steps
        m_avg = m_avg / real_MC_steps
        m_var = m_var / real_MC_steps
        susceptibility = (m_var - m_avg*m_avg)*beta(j)
        capacity = (u_var - u_avg*u_avg)*beta(j)*beta(j)

        !$omp critical
            write(unit_temperature,format_style2) KbT(j), threadID, m_avg, susceptibility, u_avg, capacity
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