
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
    implicit none

    real(pr), dimension(:), allocatable     :: Pressures, Temperatures
    real(pr), dimension(:,:), allocatable   :: positions, velocities, forces, previous_forces, Energies
!    real(pr)                               :: CPU_t_start, CPU_t_end, CPU_elapsed_time
!    character (len=:), allocatable          :: filename, prefix, file_root_word, suffix
    integer(int_huge)                       :: i, j
    integer(int_medium)                     :: unit_positions, unit_observables, unit_info

    abstract interface
        subroutine init_pos(positions)
            use precision
            implicit none
            real(pr), allocatable, intent(out) :: positions(:,:)
        end subroutine init_pos
        subroutine force_sub(positions, forces, E_potential, pressure_virial)
        use precision
            real(pr), intent(in)    :: positions(:,:)
            real(pr), intent(out)   :: forces(:,:)
            real(pr), intent(out)   :: E_potential, pressure_virial
        end subroutine
    end interface

    procedure(init_pos), pointer    :: initialize_positions => null()
    procedure(force_sub), pointer        :: get_forces => null()

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
    print*, "Structure selected: ", structure

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
            get_forces =>  get_forces_linkedList
            call create_maps(3,3,3)
            call create_links(positions)
    end select

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    if (integrator == 'velocity-Verlet') then
        print*,"-------------------- Calculating with ",integrator," --------------------"
        call get_forces(positions, forces,  Energies(1,1), pressures(1))

        open(newunit=unit_positions, file="datos/positions.xyz", status="replace")
        open(newunit=unit_observables, file="datos/observables.out", status="replace")
            if (save_transitory) call write_XYZfile(positions, velocities, 0._pr, unit_positions)
            if (save_observables) write(unit_observables,*) "## E_pot       |      E_kin      |      P       |       T"

            do i = 1, transitory_steps/rescale_steps
                do j = 1, rescale_steps
                    call update_positions_velVer(positions, velocities, forces)
                    previous_forces = forces
                    call get_forces(positions, forces, Energies(1,1), pressures(1))
                    call update_velocities_velVer(velocities, forces, previous_forces)
                    if (save_transitory) call write_XYZfile(positions, velocities, real(i,pr)*dt, unit_positions)
                end do
                call rescale_velocities(velocities, initial_Temp)
            end do

            transitory = .False.

            do i = transitory_steps + 1 , MD_steps
                call update_positions_velVer(positions, velocities, forces)
                previous_forces = forces
                call get_forces(positions, forces, Energies(1,i), pressures(i))
                call update_velocities_velVer(velocities, forces, previous_forces)
                call get_observables(velocities, Energies(2,i), pressures(i), temperatures(i))
                call write_XYZfile(positions, velocities, real(i,pr)*dt, unit_positions)
                if (save_observables) write(unit_observables, format_style) energies(:,i), pressures(i), temperatures(i)
            end do
        close(unit_positions)
        close(unit_observables)
    end if

    open(newunit=unit_info, file="datos/INFO.out", status="replace")
        write(unit_info,*) 
    close(unit_info)

end program ex1