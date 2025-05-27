
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
! Room for improvement: Name list block in a module. Declare global variables in modules
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

    real(kind=pr)                               :: pressure_virial
    real(kind=pr), dimension(:,:), allocatable  :: positions, velocities, forces, previous_forces, energies
!    real(kind=pr)                               :: CPU_t_start, CPU_t_end, CPU_elapsed_time
    character (len=:), allocatable              :: filename, prefix, file_root_word, suffix
    integer(kind=int_huge)                      :: i
    integer(int_medium)                         :: unit_positions

    abstract interface
        subroutine init_pos(positions, passed_num_atoms)
            use precision
            implicit none
            real(kind=pr), dimension(:,:), allocatable, intent(out) :: positions
            integer(kind=int_huge), intent(in)                      :: passed_num_atoms
        end subroutine init_pos
        subroutine pot(particle_distance_squared, force_contribution, E_potential, pressure_virial, potential_cutoff)
            use precision
            real(kind=pr), intent(in)               :: particle_distance_squared
            real(kind=pr), intent(out)              :: force_contribution
            real(kind=pr), intent(inout)            :: E_potential, pressure_virial
            real(kind=pr), intent(in)               :: potential_cutoff
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
    prefix = "./datos/"
    suffix = ".out"


    select case (structure)
    case ("FCC")
        initialize_positions =>  initialize_positions_FCC
    case ("BCC")
        initialize_positions =>  initialize_positions_BCC
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

    ! Define conversion factors to adimensionalize the variables
    conversion_factors = (/sigma,sigma*sqrt(molar_mass/epsilon),epsilon/Boltzmann_constant, epsilon/) ! distance, time, temperature, energy
    lattice_constant = lattice_constant/conversion_factors(1)
    start_time = start_time/conversion_factors(2)
    end_time = end_time/conversion_factors(2)
    dt = dt/conversion_factors(2)
    initial_Temp = initial_Temp/conversion_factors(3)

    periodicity = cell_dim*lattice_constant

    allocate(Energies(2,0:nint((end_time-start_time)/dt))) ! energies = (E_potential, E_kinetic)

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    if (do_velocity_verlet) then
        print*,"-------------------- Calculating with velocity-Verlet --------------------"
        call initialize_positions(positions, num_atoms)
        allocate(velocities(size(positions,1),size(positions,2)))
        allocate(forces(size(positions,1),size(positions,2)))
        call initialize_velocities(velocities, initial_Temp)
        call get_forces(positions, forces, potential,  Energies(1,0), pressure_virial, radius_cutoff)
        call get_E_kinetic(velocities,Energies(2,0))

        open(newunit=unit_positions, file="datos/positions.out", status="replace")
            do i = 1, nint((end_time-start_time)/dt)
                call update_positions_velVer(positions, velocities, forces, dt)
                previous_forces = forces
                call get_forces(positions, forces, potential, Energies(1,i), pressure_virial, radius_cutoff)
                call update_velocities_velVer(velocities, forces, previous_forces, dt)
                call get_E_kinetic(velocities, Energies(2,i))
                call write_to_file(positions, i*dt, unit_positions)
                !print*, i*dt, energies(:,i)
            end do
        close(unit_positions)
    end if

end program ex1