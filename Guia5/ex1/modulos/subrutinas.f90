MODULE subrutinas
    use precision
    use constantes
    use formats
    use funciones
    use mzranmod
    use omp_lib
    implicit none

    private     radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor
    real(pr)                :: radius_cutoff_squared, pair_corr_cutoff_sqr, potential_cutoff, Temp_factor, Pressure_factor

    integer(int_medium)     :: cell_dim(3)
    integer(int_large)      :: num_atoms, pair_corr_bins
    real(pr)                :: conversion_factors(5), periodicity(3)
    real(pr)                :: lattice_constant, dt, dtdt
    real(pr)                :: sigma, epsilon
    real(pr)                :: radius_cutoff, pair_corr_cutoff, dr, initial_Temp, density, molar_mass
    logical                 :: transitory, do_pair_correlation
    procedure(pot), pointer :: potential => null()

    abstract interface ! Intended to allow for the implementation of a different potential later on
        subroutine pot(particle_distance_squared, particle_separation, force_contribution, E_potential, pressure_virial &
            , potential_cutoff)
            use precision
            real(pr), intent(in)               :: particle_distance_squared, particle_separation(3)
            real(pr), intent(out)              :: force_contribution(3)
            real(pr), intent(inout)            :: E_potential, pressure_virial
            real(pr), intent(in)               :: potential_cutoff
        end subroutine pot
    end interface

contains

subroutine parameters_initialization()

    ! Define conversion factors to adimensionalize the variables
    conversion_factors(1) = sigma                                               ! Distance
    conversion_factors(2) = sigma*sqrt(molar_mass/(epsilon*Avogadro_number))    ! Time
    conversion_factors(3) = epsilon/Boltzmann_constant                          ! Temperature
    conversion_factors(4) = epsilon                                             ! Energy
    conversion_factors(5) = sqrt((epsilon*Avogadro_number)/molar_mass)          ! Velocity

    lattice_constant = lattice_constant/conversion_factors(1)
    initial_Temp = initial_Temp/conversion_factors(3)
    periodicity = cell_dim*lattice_constant
    radius_cutoff = radius_cutoff/conversion_factors(1)

    dt = dt/conversion_factors(2)
    dtdt = dt*dt

    radius_cutoff_squared = radius_cutoff*radius_cutoff
    potential_cutoff = Lennard_Jones_potencial(radius_cutoff_squared)

    transitory = .True.    ! Flag to avoid calculations and saving variables during the transitory steps
    Temp_factor = 2.0_pr / (3.0_pr * num_atoms)
    Pressure_factor = 1._pr / (3._pr * product(cell_dim)*lattice_constant*lattice_constant*lattice_constant)

    if(do_pair_correlation) then
        pair_corr_cutoff_sqr = pair_corr_cutoff*pair_corr_cutoff
        dr = pair_corr_cutoff/real(pair_corr_bins,pr)
    end if

end subroutine parameters_initialization

subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                :: fmt  ! Format descriptor
    character(len=12)               :: x1   ! Temporary string for formatted real
    character(len=*), intent(in)    :: prefix, suffix
    character(len=*)                :: filename
    real(kind=pr), intent(in)       :: num  ! Input real number

    fmt = '(I10)'  ! Format integer
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

subroutine initialize_positions_random(positions) ! Not debugged
    real(pr), allocatable, intent(out) :: positions(:,:)
    integer                            :: i

    do i = 1, size(cell_dim)
        num_atoms = num_atoms * cell_dim(i)
    end do

    ! Need to add a check to avoid collisions
    positions = reshape( [ (rmzran(), i = 1, size(positions)) ], shape(positions) )  ! Need to redefine this to get random positions in supercell, not unit cell

end subroutine initialize_positions_random

subroutine initialize_positions_FCC(positions)
    real(pr), allocatable, intent(out) :: positions(:,:)
    integer(int_large)                 :: supercell_atoms
    integer(int_large)                 :: atom_id
    integer(int_medium)                :: h, k, l, b
    integer(int_small)                 :: i

    ! FCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,4) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.0_pr, &  ! Atom 2
        0.5_pr, 0.0_pr, 0.5_pr, &  ! Atom 3
        0.0_pr, 0.5_pr, 0.5_pr   & ! Atom 4
    ], [3,4])

    supercell_atoms = 4
    do i = 1_int_small, size(cell_dim)
        supercell_atoms = supercell_atoms * cell_dim(i)
    end do
    if (supercell_atoms/=num_atoms .and. num_atoms/=0) then
        print'(a,I5,a,I5)', "Number of atoms mismatch: selected number of atoms = ", num_atoms &
            , " and atoms in the FCC supercell = ",supercell_atoms
        print*, "Ignoring number of atoms specified"
        num_atoms = supercell_atoms
    end if
    allocate(positions(3, num_atoms))

    atom_id = 1
    do h = 0_int_medium, cell_dim(1)-1_int_medium
        do k = 0_int_medium, cell_dim(2)-1_int_medium
            do l = 0_int_medium, cell_dim(3)-1_int_medium
                do b = 1, 4
                    positions(:, atom_id) = [real(h, pr), real(k, pr), real(l, pr)] + basis(:, b)
                    atom_id = atom_id + 1
                end do
            end do
        end do
    end do

    positions = positions*lattice_constant

end subroutine initialize_positions_FCC

subroutine initialize_positions_BCC(positions)
    real(pr), allocatable, intent(out) :: positions(:,:)
    integer(int_large)                 :: supercell_atoms
    integer(int_large)                 :: atom_id
    integer(int_medium)                :: h, k, l, b
    integer(int_small)                 :: i

    ! BCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,2) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.5_pr  &  ! Atom 2
    ], [3,2])

    supercell_atoms = 2
    do i = 1_int_small, size(cell_dim)
        supercell_atoms = supercell_atoms * cell_dim(i)
    end do
    if (supercell_atoms/=num_atoms .and. num_atoms/=0) then
        print'(a,I5,a,I5)', "Number of atoms mismatch: selected number of atoms = ", num_atoms &
            , " and atoms in the BCC supercell = ",supercell_atoms
        print*, "Ignoring number of atoms specified"
        num_atoms = supercell_atoms
    end if
    allocate(positions(3, num_atoms))

    atom_id = 1
    do h = 0_int_medium, cell_dim(1)-1_int_medium
        do k = 0_int_medium, cell_dim(2)-1_int_medium
            do l = 0_int_medium, cell_dim(3)-1_int_medium
                do b = 1, 4
                    positions(:, atom_id) = [real(h, pr), real(k, pr), real(l, pr)] + basis(:, b)
                    atom_id = atom_id + 1
                end do
            end do
        end do
    end do

end subroutine initialize_positions_BCC

subroutine initialize_velocities(velocities, initial_Temp)
    real(pr), allocatable, intent(out)          :: velocities(:,:)
    real(pr), intent(in)                        :: initial_Temp
    integer                                     :: i

    velocities = reshape( [ (rmzran() - 0.5_pr, i = 1, size(velocities)) ], shape(velocities) )

    call rescale_velocities(velocities, initial_Temp)

end subroutine initialize_velocities

subroutine rescale_velocities(velocities, initial_Temp)
    real(pr), allocatable, intent(inout)          :: velocities(:,:)
    real(pr), intent(in)                        :: initial_Temp
    real(pr), dimension(size(velocities,1))     :: velocity_average, velocity_variance, scaling_factor
    integer                                     :: i

    velocity_average = sum(velocities,2)
    velocity_variance = sum(velocities*velocities,2)

    velocity_average = velocity_average/real(size(velocities,2))
    velocity_variance = velocity_variance/real(size(velocities,2))
    scaling_factor = sqrt( 3._pr*initial_Temp/velocity_variance )

    do i = 1, size(velocities,1)
        velocities(i,:) = (velocities(i,:)-velocity_average(i))*scaling_factor(i)
    end do

end subroutine rescale_velocities

subroutine get_forces_allVSall(positions, forces, E_potential, pressure_virial, pair_corr)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer(int_huge)       :: i, j

    forces = 0._pr
    if (transitory) then; E_potential = 0.0; pressure_virial = 0.0 ; end if

    !$omp parallel private(j, i) &
    !$omp shared(positions, forces, num_atoms) &
    !$omp reduction(+:E_potential, pressure_virial, pair_corr)

        !$omp do schedule(dynamic)
        do i=1,num_atoms-1
            do j = i+1, num_atoms
                call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential, &
                    pressure_virial, pair_corr)
            end do
        end do
        !$omp end do

    !$omp end parallel


end subroutine get_forces_allVSall

subroutine get_force_contribution(particle1_position, particle2_position, particle1_forces, particle2_forces, E_potential &
    , pressure_virial, pair_corr)
    real(pr), dimension(3), intent(in)      :: particle1_position, particle2_position
    real(pr), dimension(3), intent(out)     :: particle1_forces, particle2_forces
    real(pr), intent(out)                   :: E_potential, pressure_virial
    real(pr), intent(inout)                 :: pair_corr(:)
    real(pr), dimension(3)                  :: particle_separation, force_contribution
    real(pr)                                :: particle_distance_squared

    particle_separation = particle1_position - particle2_position ! Separation vector
    particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
    particle_distance_squared = sum(particle_separation*particle_separation)
    if (particle_distance_squared <= radius_cutoff_squared) then
        call potential(particle_distance_squared, particle_separation, force_contribution, E_potential &
            , pressure_virial, potential_cutoff)
        particle1_forces = particle1_forces + force_contribution
        particle2_forces = particle2_forces - force_contribution
    endif

    if (do_pair_correlation) call update_pair_correlation(particle_distance_squared, pair_corr)

end subroutine get_force_contribution

subroutine Lennard_Jones(particle_distance_squared, particle_separation,  force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_squared, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff
    real(pr)                   :: r2inv, r6inv, force_magnitude

    r2inv = 1._pr/particle_distance_squared
    r6inv = r2inv*r2inv*r2inv
    force_magnitude = 48._pr*r2inv*r6inv*(r6inv-0.5_pr)
    force_contribution = force_magnitude*particle_separation

    E_potential = E_potential + 4._pr*r6inv*(r6inv-1._pr) - potential_cutoff
    pressure_virial = pressure_virial + particle_distance_squared*force_magnitude

end subroutine Lennard_Jones

subroutine update_positions_velVer(positions, velocities, forces)
    real(pr), dimension(:,:), intent(inout)    :: positions, velocities, forces

    positions = positions + velocities*dt + forces*dtdt*0.5_pr

    ! Apply periodic boundary conditions
    !positions = mod(positions, spread(periodicity, dim=2, ncopies=size(positions,2)))
    positions(1,:) = modulo(positions(1,:), periodicity(1))
    positions(2,:) = modulo(positions(2,:), periodicity(2))
    positions(3,:) = modulo(positions(3,:), periodicity(3))

end subroutine update_positions_velVer

subroutine update_velocities_velVer(velocities, forces, previous_forces)
    real(pr), dimension(:,:), intent(inout)    :: velocities, forces, previous_forces

    velocities = velocities + dt * 0.5_pr*(previous_forces + forces)

end subroutine update_velocities_velVer

subroutine update_pair_correlation(particle_distance_squared, pair_corr)
    real(pr), intent(in)        :: particle_distance_squared
    real(pr), intent(inout)     :: pair_corr(:)
    integer(int_large)          :: radius_bin

    if (particle_distance_squared <= pair_corr_cutoff_sqr) then
        radius_bin = int(sqrt(particle_distance_squared) / dr) + 1
        pair_corr(radius_bin) = pair_corr(radius_bin) + 1._pr
    end if

end subroutine update_pair_correlation

subroutine normalize_pair_correlation(pair_corr, MD_steps)
    real(pr), intent(inout)     :: pair_corr(:)
    integer                     :: i, MD_steps
    real(pr)                    :: r_lower, r_upper, shell_vol, ideal_pair_corr

    do i = 1, pair_corr_bins
        r_lower = real(i-1,pr) * dr
        r_upper = r_lower + dr
        shell_vol = (4.0_pr*pi / 3.0_pr) * (r_upper**3 - r_lower**3)
        ideal_pair_corr = density * shell_vol * num_atoms
        pair_corr(i) = pair_corr(i) / ideal_pair_corr / MD_steps
    end do

end subroutine normalize_pair_correlation

subroutine get_observables(velocities, E_kinetic, Pressure, Temperature)
    real(pr), intent(in)   :: velocities(:,:)
    real(pr), intent(out)  :: E_kinetic, Pressure, Temperature

    E_kinetic = 0.5_pr*sum(velocities*velocities)
    Temperature = Temp_factor*E_kinetic
    Pressure = density*Temperature + Pressure*Pressure_factor

end subroutine get_observables

subroutine get_stats(measurements, variance, stddev, average)
    real(pr), intent(in)                :: measurements(:)
    real(pr), optional, intent(out)     :: variance, stddev, average
    real(pr)                            :: avg, sqr_avg

    if (present(average) .or. present(variance) .or. present(stddev)) then
        avg = sum(measurements)/size(measurements)
        sqr_avg = sum(measurements*measurements)/size(measurements)
    end if

    if  (present(average))   average  = avg
    if  (present(variance))  variance = sqr_avg - avg*avg
    if  (present(stddev))     stddev   = sqrt(sqr_avg - avg*avg)

end subroutine get_stats


!##################################################################################################
!     Old/Unused
!##################################################################################################

subroutine get_forces_old(positions, forces, E_potential, pressure_virial)
    real(pr), dimension(:,:), intent(out)   :: positions, forces
    real(pr), intent(out)                   :: E_potential, pressure_virial
    real(pr), dimension(3)                  :: particle1_position, particle2_position, particle_separation
    real(pr)                                :: particle_distance_squared, force_contribution(3)
    integer(int_huge)                       :: i, j

    forces = 0._pr
    if (transitory) then; E_potential = 0.0; pressure_virial = 0.0 ; end if

    do i=1,size(positions,2)-1
        particle1_position = positions(:,i)
        do j = i+1,size(positions,2)
            particle2_position = positions(:,j)
            particle_separation = particle1_position - particle2_position ! Separation vector
            particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
            particle_distance_squared = sum(particle_separation*particle_separation)
            if (particle_distance_squared <= radius_cutoff_squared) then
                call potential(particle_distance_squared, particle_separation, force_contribution, E_potential &
                    , pressure_virial, potential_cutoff)
                forces(:,i) = forces(:,i) + force_contribution
                forces(:,j) = forces(:,j) - force_contribution
            endif
        end do
    end do

end subroutine get_forces_old



END MODULE