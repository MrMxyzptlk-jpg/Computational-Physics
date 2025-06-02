
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
! Room for improvement: add check to the random initialization of positions in order to avoid collisions
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
    implicit none

    real(pr), dimension(:), allocatable     :: Pressures, Temperatures
    real(pr), dimension(:,:), allocatable   :: positions, velocities, forces, previous_forces, Energies
!    real(pr)                               :: CPU_t_start, CPU_t_end, CPU_elapsed_time
!    character (len=:), allocatable          :: filename, prefix, file_root_word, suffix
    integer(int_huge)                       :: i, j, N_iterations
    integer(int_medium)                     :: unit_positions, unitnum

    abstract interface
        subroutine init_pos(positions)
            use precision
            implicit none
            real(pr), allocatable, intent(out) :: positions(:,:)
        end subroutine init_pos

        subroutine pot(particle_distance_squared, force_contribution, E_potential, pressure_virial, potential_cutoff)
            use precision
            real(pr), intent(in)               :: particle_distance_squared
            real(pr), intent(out)              :: force_contribution
            real(pr), intent(inout)            :: E_potential, pressure_virial
            real(pr), intent(in)               :: potential_cutoff
        end subroutine pot
    end interface

    procedure(init_pos), pointer    :: initialize_positions => null()
    procedure(pot), pointer         :: potential => null()

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
            if (density>0) lattice_constant = (2._pr/density)**(1._pr/3._pr)
        case ("BCC")
            initialize_positions =>  initialize_positions_BCC
            if (density>0) lattice_constant = (4._pr/density)**(1._pr/3._pr)
        case ("random")
            initialize_positions =>  initialize_positions_random
        case default
            initialize_positions =>  initialize_positions_random
    end select

    print*, "Structure selected: ", structure

    select case (type)
        case ("lannard_jones")
            potential =>  Lennard_Jones
        case default
            potential =>  Lennard_Jones
    end select

    ! Initialize random number generator
    call mzran_init()

    call parameters_initialization()

    !N_iterations = MD_steps - transitory_steps

    allocate(Energies(2,1:MD_steps)) ! energies = (E_potential, E_kinetic)
    allocate(Pressures(1:MD_steps), Temperatures(1:MD_steps))

    call initialize_XYZ_data()

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    if (integrator == 'velocity-Verlet') then
        print*,"-------------------- Calculating with ",integrator," --------------------"
        call initialize_positions(positions)
        allocate(velocities(size(positions,1),size(positions,2)))
        allocate(forces(size(positions,1),size(positions,2)))
        call initialize_velocities(velocities, initial_Temp)
        call get_forces(positions, forces, potential,  Energies(1,1), pressures(1))

        open(newunit=unit_positions, file="datos/positions.xyz", status="replace")
            call write_to_XYZfile(positions, 0._pr, unit_positions)

            do i = 1, transitory_steps/rescale_steps
                do j = 1, rescale_steps
                    call update_positions_velVer(positions, velocities, forces)
                    previous_forces = forces
                    call get_forces(positions, forces, potential, Energies(1,1), pressures(1))
                    call update_velocities_velVer(velocities, forces, previous_forces)
                    call write_to_XYZfile(positions, real(i,pr)*dt, unit_positions)
                end do
                call rescale_velocities(velocities, initial_Temp)
            end do

            transitory = .False.

            do i = transitory_steps + 1 , MD_steps
                call update_positions_velVer(positions, velocities, forces)
                previous_forces = forces
                call get_forces(positions, forces, potential, Energies(1,i), pressures(i))
                call update_velocities_velVer(velocities, forces, previous_forces)
                call get_observables(velocities, Energies(2,i), pressures(i), temperatures(i))
                call write_to_XYZfile(positions, real(i,pr)*dt, unit_positions)
            end do
        close(unit_positions)
    end if

end program ex1