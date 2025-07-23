MODULE constantsMod
    use precisionMod
    implicit none

    !real (kind=pr), parameter  :: Boltzmann_constant   = 1.380649e-23  ! In SI units
    real (pr), parameter        :: pi = acos(-1.0_pr)
    real(pr), parameter         :: twoPi = 2.0_pr*pi
    real(pr), parameter         :: twoPi_sqr = twoPi**2
    real(pr), parameter         :: fourpi = 4._pr*pi

END MODULE constantsMod