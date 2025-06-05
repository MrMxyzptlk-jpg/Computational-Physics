MODULE linkedLists
    use precision
    use funciones
    use subrutinas
    implicit none

    private  HEAD, LIST, map, N_linkedCells
    real(pr)                :: side_linkCell(3), side_inv_linkCell(3)
    integer                 :: N_linkedCells, dim_linkCell(3)
    integer, allocatable    :: HEAD(:), LIST(:), map(:) ! M: debe ser menor o igual que L/int(L/rc

contains

subroutine create_maps(dimx_linkCell, dimy_linkCell, dimz_linkCell)
    integer, intent(in) :: dimx_linkCell, dimy_linkCell, dimz_linkCell
    integer(int_large)  :: ix, iy, iz, imap

    ! Initializing global variables
    dim_linkCell      = (/dimx_linkCell, dimy_linkCell, dimz_linkCell/)
    N_linkedCells     = product(dim_linkCell)
    side_linkCell      = lattice_constant*dim_linkCell / real(dim_linkCell, pr)
    side_inv_linkCell = 1._pr / side_linkCell  ! Inverse cell size
    allocate(HEAD(N_linkedCells), LIST(num_atoms), map(13*N_linkedCells))

    do ix=1, dim_linkCell(1)
        do iy=1, dim_linkCell(2)
            do iz=1, dim_linkCell(3)
                imap = ( index_cell (ix, iy, iz, dim_linkCell) - 1 ) * 13
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
    HEAD(1:N_linkedCells) = 0

    do i = 1, size(list)
        ! Get cell indices (0-based) and then periodic boundary correction (modulo M)
        position_index = (/(mod(int( positions(j,i) * side_inv_linkCell(j) ), dim_linkCell(j)), j = 1, 3)/)

        ! Compute cell index (1-based Fortran indexing)
        jcell = 1 + position_index(1) + position_index(2)*dim_linkCell(1) + position_index(3)*dim_linkCell(1)*dim_linkCell(2)

        ! Insert particle i at the head of the list for this cell
        LIST(i) = HEAD(jcell)
        HEAD(jcell) = i
    end do
end subroutine create_links

subroutine get_forces_linkedList(positions, forces, E_potential, pressure_virial)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: forces(:,:)
    real(pr), intent(out)   :: E_potential, pressure_virial
    integer                 :: i, j, icell, jcell, jcell0, neighbor

    do icell = 1, N_linkedCells ! Go through all cells
        i = head(icell)
        do while (i /= 0)
            j = list(i)
            do while (j /= 0) ! All pairs in the cell
                call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential, pressure_virial)
                j = list(j)
            end do
            jcell0 = 13*(icell - 1)

            do neighbor = 1, 13 ! Go through all neighbor cells
                jcell = map(jcell0 + neighbor)
                j = head(jcell)

                do while (j /= 0) ! For all particles in a neighbor cell
                    call get_force_contribution(positions(:,i), positions(:,j), forces(:,i), forces(:,j), E_potential &
                        , pressure_virial)
                    j = list(j)
                end do

            end do
            i = list(i)
        end do
    end do
end subroutine get_forces_linkedList

END MODULE linkedLists