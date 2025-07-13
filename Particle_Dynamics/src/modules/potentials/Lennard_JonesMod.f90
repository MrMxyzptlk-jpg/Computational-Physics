MODULE Lennard_JonesMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    implicit none

CONTAINS

function Lennard_Jones_potential(index1, index2, particle_distance_sqr)
    real(pr), intent(in)   :: particle_distance_sqr
    integer, intent(in)    :: index1, index2
    real(pr)               :: Lennard_Jones_potential
    real(pr)               :: r2inv, r6inv

    r2inv = 1._pr/particle_distance_sqr
    r6inv = r2inv*r2inv*r2inv
    Lennard_Jones_potential =  4.0_pr *r6inv * (r6inv - 1.0_pr) - potential_cutoff

end function Lennard_Jones_potential

subroutine Lennard_Jones(index1, index2, particle_distance_sqr, particle_separation,  force_contribution, E_potential &
    , pressure_virial, potential_cutoff)
    real(pr), intent(in)    :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)   :: force_contribution(3)
    real(pr), intent(inout) :: E_potential, pressure_virial
    real(pr), intent(in)    :: potential_cutoff
    integer, intent(in)     :: index1, index2
    real(pr)                :: r2inv, r6inv, force_magnitude

    r2inv = 1._pr/particle_distance_sqr
    r6inv = r2inv*r2inv*r2inv
    force_magnitude = 48._pr*r2inv*r6inv*(r6inv-0.5_pr)
    force_contribution = force_magnitude*particle_separation

    if (measure .and. save_observables) then
        E_potential = E_potential + 4._pr*r6inv*(r6inv-1._pr) - potential_cutoff
        pressure_virial = pressure_virial + particle_distance_sqr*force_magnitude
    end if

end subroutine Lennard_Jones

END MODULE Lennard_JonesMod