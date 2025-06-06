MODULE constantes
    use precision
    implicit none

    real (kind=pr), parameter   :: Boltzmann_constant   = 1.380649e-23  ! In SI units
    real (kind=pr), parameter   :: Avogadro_number      = 6.02214076e23
    real (kind=pr), parameter   :: pi                   = acos(-1.0_pr)

END MODULE