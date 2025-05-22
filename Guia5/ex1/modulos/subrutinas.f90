MODULE subrutinas
use precision
use constantes
use formats
use funciones
use mzranmod
implicit none
    integer(int_medium)      :: cell_dim(3)
    integer(int_huge)        :: num_atoms

contains

subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                        :: fmt  ! Format descriptor
    character(len=12)                       :: x1   ! Temporary string for formatted real
    character(len=:), allocatable           :: prefix, suffix, filename
    integer(int_huge), intent(in)    :: num  ! Input real number

    fmt = '(I10)'  ! Format integer
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

subroutine write_to_file(Y, time, dx, conversion_factors, unitnum)
    real (pr), dimension (:), intent (in)  :: Y, conversion_factors
    real (pr), intent (in)                 :: dx, time
    integer (int_medium)                   :: unitnum
    integer                                     :: i

    do i = 1, size(Y)
        write(unitnum, fmt=format_style2) time*conversion_factors(2) &
        , dx*real(i-1,pr)*conversion_factors(1), Y(i)
    end do

end subroutine write_to_file

subroutine initialize_positions_FCC(positions, passed_num_atoms)
    real(pr), dimension(:,:), allocatable, intent(out) :: positions
    integer(int_huge), intent(in)                      :: passed_num_atoms
    integer(int_huge)                                  :: atom_id
    integer(int_medium)                                :: h,k,l,b
    integer(int_small)                                 :: i

    ! FCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,4) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.0_pr, &  ! Atom 2
        0.5_pr, 0.0_pr, 0.5_pr, &  ! Atom 3
        0.0_pr, 0.5_pr, 0.5_pr   & ! Atom 4
    ], [3,4])

    num_atoms = 4
    do i = 1_int_small, size(cell_dim)
        num_atoms = num_atoms * cell_dim(i)
    end do
    if (passed_num_atoms/=num_atoms .and. passed_num_atoms/=0) then
        print*, "Number of atoms mismatch: selected number of atoms = ", passed_num_atoms, "atoms in the BCC supercell = ",num_atoms
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

end subroutine initialize_positions_FCC

subroutine initialize_positions_BCC(positions, passed_num_atoms)
    real(pr), dimension(:,:), allocatable, intent(out) :: positions
    integer(int_huge), intent(in)                      :: passed_num_atoms
    integer(int_huge)                                  :: atom_id, num_atoms
    integer(int_medium)                                :: h,k,l,b
    integer(int_small)                                 :: i

    ! BCC as simple cubic with a basis:
    real(pr), parameter         :: basis(3,2) = reshape([ &
        0.0_pr, 0.0_pr, 0.0_pr, &  ! Atom 1
        0.5_pr, 0.5_pr, 0.5_pr  &  ! Atom 2
    ], [3,2])

    num_atoms = 2
    do i = 1_int_small, size(cell_dim)
        num_atoms = num_atoms * cell_dim(i)
    end do
    if (passed_num_atoms/=num_atoms .and. passed_num_atoms/=0) then
        print*, "Number of atoms mismatch: selected number of atoms = ", passed_num_atoms, "atoms in the BCC supercell = ",num_atoms
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
    real(pr), dimension(size(velocities,1))     :: velocity_average, velocity_variance, scaling_factor
    integer                                     :: i

    velocities = reshape( [ (rmzran() - 0.5_pr, i = 1, size(velocities)) ], shape(velocities) )

    velocity_average = sum(velocities,2)
    velocity_variance = sum(velocities*velocities,2)

    velocity_average = velocity_average/real(size(velocities,2))
    velocity_variance = velocity_variance/real(size(velocities,2))
    scaling_factor = sqrt( 3._pr*initial_Temp/velocity_variance )

    do i = 1, size(velocities,1)
        velocities(i,:) = (velocities(i,:)-velocity_average(i))*scaling_factor(i)
    end do

end subroutine initialize_velocities

subroutine initialize_positions_random(positions, passed_num_atoms)
    real(pr), allocatable, intent(out) :: positions(:,:)
    integer(int_huge), intent(in)      :: passed_num_atoms
    integer(int_huge)                  :: num_atoms
    integer                            :: i

    num_atoms = passed_num_atoms
    do i = 1, size(cell_dim)
        num_atoms = num_atoms * cell_dim(i)
    end do

    positions = reshape( [ (rmzran(), i = 1, size(positions)) ], shape(positions) )  ! Need to redefine this to get random positions in supercell, not unit cell

end subroutine initialize_positions_random

subroutine get_forces(positions, forces, potential, E_potential, pressure_virial, radius_cutoff, periodicity)
    real(pr), dimension(:,:), intent(out)  :: positions, forces
    real(pr), dimension(3), intent(in)     :: periodicity
    real(pr), intent(in)                   :: radius_cutoff
    real(pr), intent(out)                  :: E_potential, pressure_virial
    real(pr), dimension(3)                 :: particle1_position, particle2_position, particle_separation
    real(pr)                               :: particle_distance_squared, radius_cutoff_squared, force_contribution
    real(pr)                               :: potential_cutoff
    integer(int_huge)                      :: i, j

    interface
    subroutine potential(particle_distance_squared, force_contribution, E_potential, pressure_virial, potential_cutoff)
        use precision
        real(pr), intent(in)               :: particle_distance_squared
        real(pr), intent(out)              :: force_contribution
        real(pr), intent(inout)            :: E_potential, pressure_virial
        real(pr), intent(in)               :: potential_cutoff
    end subroutine potential
    end interface

    forces = 0._pr
    E_potential = 0.0
    radius_cutoff_squared = radius_cutoff*radius_cutoff
    potential_cutoff = Lennard_Jones_potencial(radius_cutoff_squared)

    do i=1,size(positions,2)-1
        particle1_position = positions(:,i)
        do j = i+1,size(positions,2)
            particle2_position = positions(:,j)
            particle_separation = particle1_position - particle2_position ! Separation vector
            particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
            particle_distance_squared = sum(particle_separation*particle_separation)
            if (particle_distance_squared <= radius_cutoff_squared) then
                call Lennard_Jones(particle_distance_squared, force_contribution, E_potential, pressure_virial, potential_cutoff)
                forces(:,i) = forces(:,i) + particle_separation*force_contribution
                forces(:,j) = forces(:,j) - particle_separation*force_contribution
            endif
        enddo
    enddo

end subroutine get_forces

subroutine get_E_kinetic(velocities, E_kinetic)
    real(pr), dimension(:,:), intent(in)   :: velocities
    real(pr), intent(out)                  :: E_kinetic

    E_kinetic = 0.5_pr*sum(velocities*velocities)

end subroutine get_E_kinetic

subroutine Lennard_Jones(particle_distance_squared, force_contribution, E_potential, pressure_virial, potential_cutoff)
    real(pr), intent(in)               :: particle_distance_squared
    real(pr), intent(out)              :: force_contribution
    real(pr), intent(inout)            :: E_potential, pressure_virial
    real(pr), intent(in)               :: potential_cutoff
    real(pr)                           :: r2inv, r6inv

    r2inv = 1._pr/particle_distance_squared
    r6inv = r2inv*r2inv*r2inv
    force_contribution = 48._pr*r2inv*r6inv*(r6inv-0.5_pr)

    E_potential = E_potential + 4._pr*r6inv*(r6inv-1._pr) - potential_cutoff
    pressure_virial = pressure_virial + particle_distance_squared*force_contribution

end subroutine Lennard_Jones

subroutine update_positions_velVer(positions, velocities, forces, dt, periodicity)
    real(pr), dimension(:,:), intent(inout)    :: positions, velocities, forces
    real(pr), dimension(3), intent(in)         :: periodicity
    real(pr), intent(in)                       :: dt
    real(pr)                                   :: dtdt

    dtdt = dt*dt ! Should be defined globally as well as dt
    positions = positions + velocities*dt + forces*dtdt*0.5_pr

    ! Apply periodic boundary conditions
    !positions = mod(positions, spread(periodicity, dim=2, ncopies=size(positions,2)))
    positions(1,:) = mod(positions(1,:), periodicity(1))

end subroutine update_positions_velVer

subroutine update_velocities_velVer(velocities, forces, previous_forces, dt)
    real(pr), dimension(:,:), intent(inout)    :: velocities, forces, previous_forces
    real(pr), intent(in)                       :: dt

    velocities = velocities + dt * 0.5_pr*(previous_forces + forces)

end subroutine update_velocities_velVer

subroutine create_maps()
    integer(int_large)                             :: ix, iy, iz
    integer(int_large)                             :: imap
    integer(int_large), dimension(1:13*cell_dim(1)*cell_dim(2)*cell_dim(3)) :: map

    do ix=1, cell_dim(1)
        do iy=1, cell_dim(2)
            do iz=1, cell_dim(3)
                imap = ( icell (ix, iy, iz, cell_dim) - 1 ) * 13
                map(imap+1)  = icell( ix+1,   iy  ,  iz   , cell_dim)
                map(imap+2)  = icell(ix +1,  iy+1 ,  iz   , cell_dim)
                map(imap+3)  = icell(  ix ,  iy+1 ,  iz   , cell_dim)
                map(imap+4)  = icell( ix-1,  iy+1 ,  iz   , cell_dim)
                map(imap+5)  = icell( ix+1,   iy  , iz-1  , cell_dim)
                map(imap+6)  = icell( ix+1,  iy+1 , iz-1  , cell_dim)
                map(imap+7)  = icell(  ix ,  iy+1 , iz-1  , cell_dim)
                map(imap+8)  = icell( ix-1,  iy+1 , iz-1  , cell_dim)
                map(imap+9)  = icell( ix+1,   iy  , iz+1  , cell_dim)
                map(imap+10) = icell( ix+1,  iy+1 , iz+1  , cell_dim)
                map(imap+11) = icell(  ix ,  iy+1 , iz+1  , cell_dim)
                map(imap+12) = icell( ix-1,  iy+1 , iz+1  , cell_dim)
                map(imap+13) = icell(  ix ,   iy  , iz+1  , cell_dim)
            enddo
        enddo
    enddo
end subroutine create_maps

subroutine create_links(positions, HEAD, LIST, side_length)
    real(pr), intent(in)                                :: side_length
    real(pr), intent(in)                                :: positions(:,:)
    integer, dimension(product(cell_dim)), intent(out)  :: head
    integer, dimension(size(positions,2)), intent(out)  :: list
    integer                                             :: i, icell
    integer(int_medium), dimension(3)                   :: position_index
    real(pr), dimension(3)                              :: cell_length_inv

    ! Initialize
    HEAD = 0
    cell_length_inv = real(cell_dim, pr)/side_length  ! Inverse cell size

    do i = 1, size(list)
        ! Get cell indices (0-based) and then periodic boundary correction (modulo M)
        position_index = mod(int( (positions(:,i) + 0.5_pr * side_length) * cell_length_inv ), cell_dim)

        ! Compute cell index (1-based Fortran indexing)
        icell = 1 + position_index(1) + position_index(2)*cell_dim(1) + position_index(3)*cell_dim(1)*cell_dim(2)

        ! Insert particle i at the head of the list for this cell
        LIST(i) = HEAD(icell)
        HEAD(icell) = i
    end do
end subroutine create_links


END MODULE