! NOT TESTED OR FULLY IMPLEMENTED
MODULE Coulomb_ReactionFieldMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    implicit none

CONTAINS

function reaction_field_potential(particle_distance_sqr) ! Coulomb potential with reaction field approximation
    real(pr), intent(in)   :: particle_distance_sqr
    real(pr)               :: reaction_field_potential

    reaction_field_potential =  1._pr/sqrt(particle_distance_sqr)

end function reaction_field_potential

subroutine reaction_field(particle_distance_sqr, particle_separation, force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff

    real(pr) :: particle_distance, particle_distance_cubed
    real(pr) :: radius_cutoff, radius_cutoff_cubed
    real(pr) :: q_i, q_j, prefac, V_coulomb, V_rf
    real(pr) :: F_coulomb, F_rf, force_magnitude
    real(pr) :: dielectric_factor

    particle_distance = sqrt(particle_distance_sqr)
    particle_distance_cubed = particle_distance * particle_distance_sqr
    radius_cutoff = sqrt(radius_cutoff_squared)
    radius_cutoff_cubed = radius_cutoff*radius_cutoff_squared

    ! Charges (unit charges)
    q_i = 1._pr
    q_j = 1._pr

    ! Prefactor: 1 / (4*pi*delta_0) in atomic units = 1
    prefac = q_i * q_j

    ! Dielectric correction factor (ε → ∞ gives 1)
    dielectric_factor = 1.0_pr

    ! Coulomb energy and force
    V_coulomb = prefac / particle_distance
    F_coulomb = prefac / particle_distance_cubed

    ! Reaction field energy and force correction
    V_rf = prefac * dielectric_factor * (particle_distance_sqr - 3._pr*radius_cutoff_squared) / (2._pr * radius_cutoff_cubed)
    F_rf = prefac * dielectric_factor * particle_distance / radius_cutoff_cubed

    ! Total force magnitude
    force_magnitude = F_coulomb + F_rf
    force_contribution = force_magnitude * particle_separation

    if (measure .and. save_observables) then
        E_potential = E_potential + V_coulomb + V_rf - potential_cutoff
        pressure_virial = pressure_virial + particle_distance_sqr * force_magnitude
    end if

end subroutine reaction_field

END MODULE Coulomb_ReactionFieldMod