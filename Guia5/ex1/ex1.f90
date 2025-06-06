
!***************************************************************
! Program: ex1.f90
! Purpose:
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
! Room for improvement: add check to the random initialization of positions in order to avoid collisions. Check cutoff_radius units (and relation to sigma).
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precision
    use funciones
    use subrutinas
    use constantes
    use mzranmod
    use parsing
    use linkedLists
    use omp_lib
    use writing2files
    implicit none

    real(pr), dimension(:), allocatable     :: Pressures, Temperatures, pair_corr
    real(pr), dimension(:,:), allocatable   :: positions, velocities, forces, previous_forces, Energies
    real(pr)                                :: CPU_t_start, CPU_t_end, CPU_elapsed_time
!    character (len=:), allocatable          :: filename, prefix, file_root_word, suffix
    integer(int_huge)                       :: i, j
    integer(int_medium)                     :: unit_positions, unit_observables
    logical                                 :: do_linkCell

    abstract interface
        subroutine init_pos(positions)
            use precision
            implicit none
            real(pr), allocatable, intent(out) :: positions(:,:)
        end subroutine init_pos
        subroutine force_sub(positions, forces, E_potential, pressure_virial, pair_corr)
        use precision
            real(pr), intent(in)    :: positions(:,:)
            real(pr), intent(out)   :: forces(:,:)
            real(pr), intent(out)   :: E_potential, pressure_virial
            real(pr), intent(inout) :: pair_corr(:)
        end subroutine
    end interface

    procedure(init_pos), pointer    :: initialize_positions => null()
    procedure(force_sub), pointer   :: get_forces => null()

!###################################################################################################
!   Set default values and parse input file
!###################################################################################################

    call set_defaults()
    call parse_input()

!##################################################################################################
!      Necessary definitions, pointers, initializations and conversion factors
!##################################################################################################

    ! Set files' name defaults
    !prefix = "./datos/"
    !suffix = ".out"


    select case (structure)
        case ("FCC")
            initialize_positions =>  initialize_positions_FCC
            if (density>0._pr) lattice_constant = (2._pr/density)**(1._pr/3._pr)
        case ("BCC")
            initialize_positions =>  initialize_positions_BCC
            if (density>0._pr) lattice_constant = (4._pr/density)**(1._pr/3._pr)
        case ("random")
            call mzran_init() ! Initialize random number generator
            initialize_positions =>  initialize_positions_random
        case default
            initialize_positions =>  initialize_positions_FCC
            if (density>0._pr) lattice_constant = (2._pr/density)**(1._pr/3._pr)
    end select

    select case (type)
        case ("lannard_jones")
            potential =>  Lennard_Jones
        case default
            potential =>  Lennard_Jones
    end select

    !N_iterations = MD_steps - transitory_steps

    call parameters_initialization()
    call initialize_XYZ_data()

    call initialize_positions(positions)
    allocate(Energies(2,1:MD_steps)) ! energies = (E_potential, E_kinetic)
    allocate(Pressures(1:MD_steps), Temperatures(1:MD_steps))
    allocate(velocities(size(positions,1),size(positions,2)))
    allocate(forces(size(positions,1),size(positions,2)))
    call initialize_velocities(velocities, initial_Temp)

    select case (summation)
        case ("all-vs-all")
            get_forces =>  get_forces_allVSall
        case ("linked-lists")
            call check_linkCell(do_linkCell)
            if (do_linkCell) then  ! Should be parsed from input
                get_forces =>  get_forces_linkedList
                call create_maps()
                call create_links(positions)
            else
                print'(a,3I3,a,3I3,a)', "Number of linked cells in each directions = (", dim_linkCell,") > L/int(L/rcut) = (" &
                    ,int(periodicity/int(periodicity/radius_cutoff)),")   --->   Using 'all-vs-all' integrator instead"
                summation = "all-vs-all"
                get_forces =>  get_forces_allVSall
            end if
    end select

    if (do_pair_correlation) then
        allocate(pair_corr(pair_corr_bins))
        pair_corr = 0._pr
    end if

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    CPU_t_start = omp_get_wtime()

    if (integrator == 'velocity-Verlet') then
        call get_forces(positions, forces,  Energies(1,1), pressures(1), pair_corr)

        open(newunit=unit_positions, file="datos/positions.xyz", status="replace")
        open(newunit=unit_observables, file="datos/observables.out", status="replace")
            if (save_transitory) call write_XYZfile(positions, velocities, 0._pr, unit_positions)
            if (save_observables) write(unit_observables,*) "## E_pot       |      E_kin      |      P       |       T"

            do i = 1, transitory_steps/rescale_steps
                do j = 1, rescale_steps
                    call update_positions_velVer(positions, velocities, forces)
                    previous_forces = forces
                    call get_forces(positions, forces, Energies(1,1), pressures(1), pair_corr)
                    call update_velocities_velVer(velocities, forces, previous_forces)
                    if (save_transitory) call write_XYZfile(positions, velocities, real(i,pr)*dt, unit_positions)
                end do
                call rescale_velocities(velocities, initial_Temp)
            end do

            transitory = .False.

            do i = 1 , MD_steps
                call update_positions_velVer(positions, velocities, forces)
                previous_forces = forces
                call get_forces(positions, forces, Energies(1,i), pressures(i), pair_corr)
                call update_velocities_velVer(velocities, forces, previous_forces)
                call get_observables(velocities, Energies(2,i), pressures(i), temperatures(i))
                call write_XYZfile(positions, velocities, real(i,pr)*dt, unit_positions)
                if (save_observables) write(unit_observables, format_style0) energies(:,i), pressures(i), temperatures(i)
            end do
        close(unit_positions)
        close(unit_observables)
    end if

    if (do_pair_correlation) then
        call normalize_pair_correlation(pair_corr, MD_steps)
        call write_pair_corr(pair_corr)
    end if

    CPU_t_end = omp_get_wtime()
    CPU_elapsed_time = CPU_t_end - CPU_t_start

    call write_output(CPU_elapsed_time)

end program ex1