MODULE funciones
    USE precision
    USE constantes
    implicit none

    contains

    !funciones de la guía 5

    function Lennard_Jones_potencial(particle_distance_squared) ! Lennard-Jones potential
        real(pr), intent(in)   :: particle_distance_squared
        real(pr)               :: Lennard_Jones_potencial
        real(pr)               :: r2inv, r6inv

        r2inv = 1._pr/particle_distance_squared
        r6inv = r2inv*r2inv*r2inv
        Lennard_Jones_potencial =  4.0_pr *r6inv * (r6inv - 1.0_pr)

    end function Lennard_Jones_potencial

    integer(int_large) function index_cell(ix,iy,iz, cell_dim)
        integer(int_medium), intent(in)    :: cell_dim(3)
        integer(int_large), intent(in)     :: ix, iy, iz
        integer(int_medium)                :: Mx, My, Mz

        Mx = cell_dim(1)
        My = cell_dim(2)
        Mz = cell_dim(3)

        index_cell = 1 + mod(ix - 1 + Mx, Mx) + mod(iy - 1 + My, My)*Mx  + mod(iz - 1 + Mz, Mz)*Mx*My

    end function index_cell

    END module