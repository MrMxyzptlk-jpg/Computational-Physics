MODULE linkedlists
    use precision
    use funciones
    use subrutinas
    use parsing, only : dim_linkCell
    use omp_lib
    implicit none

    private  head, list, map, N_linkedCells, N_neighbors
    integer, parameter      :: N_neighbors = 13
    real(pr)                :: side_linkCell(3), side_inv_linkCell(3)
    integer                 :: N_linkedCells
    integer, allocatable    :: head(:), list(:), map(:) ! M: debe ser menor o igual que L/int(L/rc

contains

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

subroutine get_forces_linkedlist(positions, forces, E_potential, pressure_virial, pair_corr)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer                 :: i, j, icell, jcell, jcell0, neighbor

    forces = 0._pr

    !$omp parallel do private(neighbor, j, jcell, jcell0, i, icell) &
    !$omp shared(positions, forces, head, map, list, N_linkedCells) &
    !$omp schedule(dynamic) reduction(+:E_potential, pressure_virial, pair_corr)
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

subroutine check_linkCell(do_linkCell)
    logical, intent(out)    :: do_linkCell
    real(pr)                :: max_cells(3)

    do_linkCell = .False.
    max_cells = periodicity/int(periodicity/radius_cutoff)

    if ((dim_linkCell(1)<=max_cells(1)).and.(dim_linkCell(2)<=max_cells(2)).and.(dim_linkCell(3)<=max_cells(3))) do_linkCell=.True.

 end subroutine check_linkCell

END MODULE linkedlists