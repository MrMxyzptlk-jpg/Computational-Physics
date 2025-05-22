MODULE funciones
    USE precision
    USE constantes
    implicit none

    contains

    !funciones de la guía 5

    function Lennard_Jones_potencial(particle_distance_squared) ! Lennard-Jones potential
        real(kind=pr), intent(in)   :: particle_distance_squared
        real(kind=pr)               :: Lennard_Jones_potencial
        real(kind=pr)               :: r2inv, r6inv

        r2inv = 1._pr/particle_distance_squared
        r6inv = r2inv*r2inv*r2inv
        Lennard_Jones_potencial =  4.0_pr *r6inv * (r6inv - 1.0_pr)

    end function Lennard_Jones_potencial

    function icell(ix,iy,iz, cell_dim)
        integer(kind=int_medium), dimension(3), intent(in)  :: cell_dim
        integer(kind=int_large), intent(in)                 :: ix, iy, iz
        integer(kind=int_medium)                            :: Mx, My, Mz
        integer(kind=int_large)                             :: icell

        Mx = cell_dim(1)
        My = cell_dim(2)
        Mz = cell_dim(3)

        icell = 1 + mod(ix - 1 + Mx, Mx) + mod(iy - 1 + My, My)*Mx  + mod(iz - 1 + Mz, Mz)*Mx*My

    end function icell

    END module