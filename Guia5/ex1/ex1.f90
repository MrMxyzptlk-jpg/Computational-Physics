
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
! Room for improvement:
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precision
    use funciones
    use subrutinas
    use constantes
    use mzranmod
    implicit none

    real(kind=pr)                               :: start_time, end_time, dt
    real(kind=pr)                               :: radius_cutoff, initial_Temp, density, lattice_constant, molar_mass
    real(kind=pr)                               :: sigma, epsilon, pressure_virial
    real(kind=pr), dimension(:,:), allocatable  :: positions, velocities, forces, previous_forces, energies
    integer(kind=int_medium), dimension(3)      :: cell_dim
    real(kind=pr), dimension(3)                 :: periodicity
    integer(kind=int_huge)                      :: passed_num_atoms
    real(kind=pr), dimension(:), allocatable    :: conversion_factors
!    real(kind=pr)                               :: CPU_t_start, CPU_t_end, CPU_elapsed_time
    character (len=:), allocatable              :: filename, prefix, file_root_word, suffix
    character (len=6)                           :: structure
    character (len=15)                          :: type
    integer(kind=int_huge)                      :: epochs, num_atoms, i
    integer                                     :: unitnum
    logical                                     :: do_velocity_verlet

    abstract interface
        subroutine init_pos(positions, passed_num_atoms, cell_dim)
            use precision
            implicit none
            real(kind=pr), dimension(:,:), allocatable, intent(out) :: positions
            integer(kind=int_medium), dimension(:), intent(in)      :: cell_dim
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

!##################################################################################################
!       Input settings
!##################################################################################################

    ! Namelist blocks
    namelist /physical/ structure, lattice_constant, density, initial_Temp, num_atoms, molar_mass, cell_dim
    namelist /calculation/ start_time, end_time, dt, radius_cutoff, epochs
    namelist /tasks/ do_velocity_verlet
    namelist /approximation/ type, sigma, epsilon

    ! DEFAULT SETTINGS
        ! Physical problems' characteristics
        structure           = "random"
        lattice_constant    = 1._pr
        cell_dim            = (/1,1,1/)
        num_atoms           = 0
        initial_Temp        = 1.1_pr
        density             = 0.8_pr
        molar_mass          = 1._pr

        !Calculation settings
        start_time      = 0._pr
        end_time        = 1800._pr
        dt              = 0.3_pr
        radius_cutoff   = 2.5_pr
        epochs          = 300

        ! Tasks
        do_velocity_verlet  = .True.

        ! Potential parameters
        sigma   = 1._pr
        epsilon = 1._pr

    ! Read from input file
    open(newunit=unitnum, file="input.nml", status="old", action="read")
        read(unitnum, nml=physical)
        read(unitnum, nml=calculation)
        read(unitnum, nml=tasks)
        read(unitnum, nml=approximation)
    close(unitnum)

!##################################################################################################
!      Necesarry definitions, pointers, initializations and conversion factors
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

    allocate(energies(2,0:nint((end_time-start_time)/dt))) ! energies = (E_potential, E_kinetic)

!##################################################################################################
!      Start of the calculations
!##################################################################################################


    if (do_velocity_verlet) then
        print*,"-------------------- Calculating with velocity-Verlet --------------------"
        call initialize_positions(positions, passed_num_atoms, cell_dim)
        allocate(velocities(size(positions,1),size(positions,2)))
        allocate(forces(size(positions,1),size(positions,2)))
        call initialize_velocities(velocities, initial_Temp)
        call get_forces(positions, forces, potential,  Energies(1,0), pressure_virial, radius_cutoff, periodicity)
        call get_E_kinetic(velocities,Energies(2,0))

        do i = 1, nint((end_time-start_time)/dt)
            call update_positions_velVer(positions, velocities, forces, dt, periodicity)
            previous_forces = forces
            call get_forces(positions, forces, potential, Energies(1,i), pressure_virial, radius_cutoff, periodicity)
            call update_velocities_velVer(velocities, forces, previous_forces, dt)
            call get_E_kinetic(velocities,Energies(2,i))
            !print*, i*dt, energies(:,i)
        end do
    end if

end program ex1