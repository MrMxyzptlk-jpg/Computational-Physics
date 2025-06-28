MODULE forcesMod
    use precisionMod
    use constantsMod
    use subroutinesMod
    use parsingMod, only : dim_linkCell, integrator
    use omp_lib
    use parametersMod
    use observablesMod
    use potentialsMod
    implicit none

    private  head, list, map, N_linkedCells, N_neighbors
    integer, parameter      :: N_neighbors = 13
    real(pr)                :: side_linkCell(3), side_inv_linkCell(3)
    integer                 :: N_linkedCells
    integer, allocatable    :: head(:), list(:), map(:) ! M: debe ser menor o igual que L/int(L/rc

CONTAINS

subroutine get_forces_allVSall(positions, forces, E_potential, pressure_virial, pair_corr)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer(int_huge)       :: i, j

    forces = 0._pr
    if (measure .and. save_observables) then; E_potential = 0.0; pressure_virial = 0.0 ; end if

    !$omp parallel private(j, i) &
    !$omp shared(positions, num_atoms) &
    !$omp reduction(+: forces, E_potential, pressure_virial, pair_corr)

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

    if (do_pair_correlation .and. .not. transitory) call update_pair_correlation(particle_distance_squared, pair_corr)

end subroutine get_force_contribution

subroutine check_linkCell(do_linkCell)
    logical, intent(out)    :: do_linkCell
    integer                 :: i, max_cells(3)

    do_linkCell = .False.

    max_cells = int(periodicity/radius_cutoff)  ! There must be less cells than the number of radius_cutoff that fit in any given direction


    if (any((/(max_cells(i) == 0 , i=1,3)/))) then
        print'(a,3I3,a)', "Radius cutoff greater than the super-cell dimensions --->   Using 'all-vs-all' integrator instead"
        return
    end if

    if (any((/(dim_linkCell(i) < 3 , i=1,3)/))) then
        print'(a,3I3,a)', "Number of linked cells in each directions = (", dim_linkCell,") < (3 3 3)   ---> "// &
            "  Using 'all-vs-all' integrator instead"
        return
    end if

    if (all(dim_linkCell <= max_cells)) then
            do_linkCell=.True.
    else
        print'(a,3I3,a,3I3,a)', "Number of linked cells in each directions = (", dim_linkCell,") > L/int(L/rcut) = (" &
            ,max_cells,")   --->   Using 'all-vs-all' integrator instead"
    end if

end subroutine check_linkCell

subroutine create_maps()
    integer(int_large)  :: ix, iy, iz, imap

    ! Initializing global variables
    N_linkedCells     = product(dim_linkCell)
    side_linkCell      = lattice_constant*dim_linkCell / real(dim_linkCell, pr)
    side_inv_linkCell = 1._pr / side_linkCell  ! Inverse cell size
    allocate(head(N_linkedCells), list(num_atoms), map(N_neighbors*N_linkedCells))

    do ix=1, dim_linkCell(1)
        do iy=1, dim_linkCell(2)
            do iz=1, dim_linkCell(3)
                imap = ( index_cell (ix, iy, iz, dim_linkCell) - 1 ) * N_neighbors
                map(imap+1)  = index_cell( ix+1,   iy  ,  iz   , dim_linkCell)
                map(imap+2)  = index_cell(ix +1,  iy+1 ,  iz   , dim_linkCell)
                map(imap+3)  = index_cell(  ix ,  iy+1 ,  iz   , dim_linkCell)
                map(imap+4)  = index_cell( ix-1,  iy+1 ,  iz   , dim_linkCell)
                map(imap+5)  = index_cell( ix+1,   iy  , iz-1  , dim_linkCell)
                map(imap+6)  = index_cell( ix+1,  iy+1 , iz-1  , dim_linkCell)
                map(imap+7)  = index_cell(  ix ,  iy+1 , iz-1  , dim_linkCell)
                map(imap+8)  = index_cell( ix-1,  iy+1 , iz-1  , dim_linkCell)
                map(imap+9)  = index_cell( ix+1,   iy  , iz+1  , dim_linkCell)
                map(imap+10) = index_cell( ix+1,  iy+1 , iz+1  , dim_linkCell)
                map(imap+11) = index_cell(  ix ,  iy+1 , iz+1  , dim_linkCell)
                map(imap+12) = index_cell( ix-1,  iy+1 , iz+1  , dim_linkCell)
                map(imap+13) = index_cell(  ix ,   iy  , iz+1  , dim_linkCell)
            end do
        end do
    end do

end subroutine create_maps

subroutine create_links(positions)
    real(pr), intent(in)    :: positions(:,:)
    integer                 :: position_index(3),  i, j, jcell

    ! Initialize
    head(1:N_linkedCells) = 0

    do i = 1, size(list)
        ! Get cell indices (0-based) and then periodic boundary correction (modulo M)
        position_index = (/(mod(int( positions(j,i) * side_inv_linkCell(j) ), dim_linkCell(j)), j = 1, 3)/)

        ! Compute cell index (1-based Fortran indexing)
        jcell = 1 + position_index(1) + position_index(2)*dim_linkCell(1) + position_index(3)*dim_linkCell(1)*dim_linkCell(2)

        ! Insert particle i at the head of the list for this cell
        list(i) = head(jcell)
        head(jcell) = i
    end do

end subroutine create_links

integer(int_large) function index_cell(ix,iy,iz, dim_linkCell)  ! For indexing Linked-Lists
    integer(int_large), intent(in)   :: dim_linkCell(3)
    integer(int_large), intent(in)   :: ix, iy, iz
    integer(int_large)               :: Mx, My, Mz

    Mx = dim_linkCell(1)
    My = dim_linkCell(2)
    Mz = dim_linkCell(3)

    index_cell = 1 + mod(ix - 1 + Mx, Mx) + mod(iy - 1 + My, My)*Mx  + mod(iz - 1 + Mz, Mz)*Mx*My

end function index_cell

subroutine get_forces_linkedlist(positions, forces, E_potential, pressure_virial, pair_corr) ! Performs worse than serial version
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer                 :: i, j, icell, jcell, jcell0, neighbor

    forces = 0._pr
    if (measure .and. save_observables) then; E_potential = 0.0; pressure_virial = 0.0 ; end if

    !$omp parallel do private(neighbor, j, jcell, jcell0, i, icell) &
    !$omp shared(positions, head, map, list, N_linkedCells) &
    !$omp schedule(dynamic) reduction(+: forces, E_potential, pressure_virial, pair_corr)
    do icell = 1, N_linkedCells ! Go through all cells
        i = head(icell)
        do while (i /= 0)
            j = list(i)
            do while (j /= 0) ! All pairs in the cell
                call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential &
                    , pressure_virial, pair_corr)
                j = list(j)
            end do
            jcell0 = N_neighbors*(icell - 1)


            do neighbor = 1, N_neighbors ! Go through all neighbor cells
                jcell = map(jcell0 + neighbor)
                j = head(jcell)

                do while (j /= 0) ! For all particles in a neighbor cell
                    call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential &
                        , pressure_virial, pair_corr)
                    j = list(j)
                end do
            end do
            i = list(i)
        end do
    end do
    !$omp end parallel do

end subroutine get_forces_linkedlist

subroutine get_forces_Ewald(positions, forces, E_potential, pressure_virial, pair_corr)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer(int_huge)       :: i, j
    real(pr)                :: box_dipole(3), surface_potential, pressure_reciprocal
    real(pr), dimension(3,num_atoms)     :: force_real, force_reciprocal
    real(pr)                   :: potential_real, potential_reciprocal

    force_reciprocal = 0._pr
    potential_reciprocal = 0._pr

    forces = 0._pr
    if (measure .and. save_observables) then; E_potential = 0.0; pressure_virial = 0.0 ; pressure_reciprocal = 0._pr; end if

    !$omp parallel private(j, i) &
    !$omp shared(positions, num_atoms) &
    !$omp reduction(+: forces, E_potential, pressure_virial, pair_corr)

        !$omp do schedule(dynamic)
        do i=1,num_atoms-1
            do j = i+1, num_atoms
                call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential, &
                    pressure_virial, pair_corr)
            end do
        end do
        !$omp end do

    !$omp end parallel

    call Coulomb_Ewald_reciprocalSpace(positions, force_reciprocal, potential_reciprocal)
    forces = forces + force_reciprocal

    ! Correction if there a net charge (non-homogenous cancelling background charge)
    if(.False.) then ! NOT DEBUGGED
        forces = forces - pi/(0.75_pr*volume) * spread(sum(positions,2)-0.5_pr*real(num_atoms,pr)*periodicity,2,num_atoms)
        if (measure .and. save_observables) then
            box_dipole = sum(positions,2) ! Net box dipole, when all charges are equal. Else the charges q must multiply each position
            surface_potential  = (pi/(1.5_pr*volume)) * sum(box_dipole*box_dipole)   ! Surface term
            E_potential = E_potential + surface_potential
        end if
    end if

    if (measure .and. save_observables) then
        E_potential = E_potential + potential_reciprocal
        pressure_virial = E_potential
    end if

end subroutine get_forces_Ewald

!##################################################################################################
!     Not used / Not implemented
!##################################################################################################

!subroutine get_forces_linkedlist_serial(positions, forces, E_potential, pressure_virial, pair_corr)
!    real(pr), intent(in)    :: positions(:,:)
!    real(pr), intent(out)   :: forces(:,:)
!    real(pr), intent(out)   :: E_potential, pressure_virial
!    real(pr), intent(inout) :: pair_corr(:)
!    integer                 :: i, j, icell, jcell, jcell0, neighbor
!
!    forces = 0._pr
!    if (measure) then; E_potential = 0.0; pressure_virial = 0.0 ; end if
!
!    do icell = 1, N_linkedCells ! Go through all cells
!        i = head(icell)
!        do while (i /= 0)
!            j = list(i)
!            do while (j /= 0) ! All pairs in the cell
!                call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential &
!                    , pressure_virial, pair_corr)
!                j = list(j)
!            end do
!            jcell0 = N_neighbors*(icell - 1)
!
!
!            do neighbor = 1, N_neighbors ! Go through all neighbor cells
!                jcell = map(jcell0 + neighbor)
!                j = head(jcell)
!
!                do while (j /= 0) ! For all particles in a neighbor cell
!                    call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential &
!                        , pressure_virial, pair_corr)
!                    j = list(j)
!                end do
!            end do
!            i = list(i)
!        end do
!    end do
!
!end subroutine get_forces_linkedlist_serial

!subroutine get_pair_correlation_linkedlist(positions, pair_corr)
!    real(pr), intent(in)    :: positions(:,:)
!    real(pr), intent(inout) :: pair_corr(:)
!    real(pr)                :: particle_distance_squared
!    integer                 :: i, j, icell, jcell, jcell0, neighbor
!
!    if (measure) then
!        !$omp parallel do private(neighbor, j, jcell, jcell0, i, icell, particle_distance_squared) &
!        !$omp shared(positions, head, map, list, N_linkedCells) &
!        !$omp schedule(dynamic) reduction(+: pair_corr)
!        do icell = 1, N_linkedCells ! Go through all cells
!            i = head(icell)
!            do while (i /= 0)
!                j = list(i)
!                do while (j /= 0) ! All pairs in the cell
!                    call get_distance_squared(positions(:,i), positions(:,j), particle_distance_squared)
!                    call update_pair_correlation(particle_distance_squared, pair_corr)
!                    j = list(j)
!                end do
!                jcell0 = N_neighbors*(icell - 1)
!
!
!                do neighbor = 1, N_neighbors ! Go through all neighbor cells
!                    jcell = map(jcell0 + neighbor)
!                    j = head(jcell)
!
!                    do while (j /= 0) ! For all particles in a neighbor cell
!                        call get_distance_squared(positions(:,i), positions(:,j), particle_distance_squared)
!                        call update_pair_correlation(particle_distance_squared, pair_corr)
!                        j = list(j)
!                    end do
!                end do
!                i = list(i)
!            end do
!        end do
!        !$omp end parallel do
!    end if
!
!end subroutine get_pair_correlation_linkedlist

!subroutine get_forces_allVSall_serial(positions, forces, E_potential, pressure_virial, pair_corr) ! Worst performing algorithm
!    real(pr), intent(in)    :: positions(:,:)
!    real(pr), intent(out)   :: forces(:,:)
!    real(pr), intent(out)   :: E_potential, pressure_virial
!    real(pr), intent(inout) :: pair_corr(:)
!    integer(int_huge)       :: i, j
!
!    forces = 0._pr
!    if (measure) then; E_potential = 0.0; pressure_virial = 0.0 ; end if
!
!    do i=1,num_atoms-1
!        do j = i+1, num_atoms
!            call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential, &
!                pressure_virial, pair_corr)
!        end do
!    end do
!
!end subroutine get_forces_allVSall_serial

END MODULE forcesMod