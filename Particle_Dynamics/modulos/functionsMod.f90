MODULE functionsMod
    USE precisionMod
    USE constantsMod
    implicit none

    procedure(pot_func), pointer :: potential_function => null()

    abstract interface
        function pot_func(particle_distance_squared)
            use precisionMod
            real(pr), intent(in)    :: particle_distance_squared
            real(pr)                :: pot_func
        end function pot_func
    end interface
    contains

integer(int_large) function index_cell(ix,iy,iz, dim_linkCell)
    integer(int_large), intent(in)   :: dim_linkCell(3)
    integer(int_large), intent(in)   :: ix, iy, iz
    integer(int_large)               :: Mx, My, Mz

    Mx = dim_linkCell(1)
    My = dim_linkCell(2)
    Mz = dim_linkCell(3)

    index_cell = 1 + mod(ix - 1 + Mx, Mx) + mod(iy - 1 + My, My)*Mx  + mod(iz - 1 + Mz, Mz)*Mx*My

end function index_cell

END module functionsMod