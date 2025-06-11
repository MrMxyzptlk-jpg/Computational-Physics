MODULE funciones
    USE precision
    USE constantes
    implicit none

    contains

    !funciones de la guía 5

    function Lennard_Jones_potential(particle_distance_squared) ! Lennard-Jones potential
        real(pr), intent(in)   :: particle_distance_squared
        real(pr)               :: Lennard_Jones_potential
        real(pr)               :: r2inv, r6inv

        r2inv = 1._pr/particle_distance_squared
        r6inv = r2inv*r2inv*r2inv
        Lennard_Jones_potential =  4.0_pr *r6inv * (r6inv - 1.0_pr)

    end function Lennard_Jones_potential

    integer(int_large) function index_cell(ix,iy,iz, dim_linkCell)
        integer(int_large), intent(in)   :: dim_linkCell(3)
        integer(int_large), intent(in)   :: ix, iy, iz
        integer(int_large)               :: Mx, My, Mz

        Mx = dim_linkCell(1)
        My = dim_linkCell(2)
        Mz = dim_linkCell(3)

        index_cell = 1 + mod(ix - 1 + Mx, Mx) + mod(iy - 1 + My, My)*Mx  + mod(iz - 1 + Mz, Mz)*Mx*My

    end function index_cell

    END module