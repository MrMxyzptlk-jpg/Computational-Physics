MODULE potentialsMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    implicit none

    procedure(pot), pointer :: potential => null()

    abstract interface ! Intended to allow for the implementation of a different potential later on
        subroutine pot(particle_distance_squared, particle_separation, force_contribution, E_potential, pressure_virial &
            , potential_cutoff)
            use precisionMod
            real(pr), intent(in)               :: particle_distance_squared, particle_separation(3)
            real(pr), intent(out)              :: force_contribution(3)
            real(pr), intent(inout)            :: E_potential, pressure_virial
            real(pr), intent(in)               :: potential_cutoff
        end subroutine pot
    end interface

CONTAINS

function Lennard_Jones_potential(particle_distance_squared) ! Lennard-Jones potential
    real(pr), intent(in)   :: particle_distance_squared
    real(pr)               :: Lennard_Jones_potential
    real(pr)               :: r2inv, r6inv

    r2inv = 1._pr/particle_distance_squared
    r6inv = r2inv*r2inv*r2inv
    Lennard_Jones_potential =  4.0_pr *r6inv * (r6inv - 1.0_pr)

end function Lennard_Jones_potential

function reaction_field_potential(particle_distance_squared) ! Coulomb potential with reaction field approximation (charges should be initialized outside)
    real(pr), intent(in)   :: particle_distance_squared
    real(pr)               :: reaction_field_potential

    reaction_field_potential =  1._pr/sqrt(particle_distance_squared)

end function reaction_field_potential

subroutine Lennard_Jones(particle_distance_squared, particle_separation,  force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_squared, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff
    real(pr)                   :: r2inv, r6inv, force_magnitude

    r2inv = 1._pr/particle_distance_squared
    r6inv = r2inv*r2inv*r2inv
    force_magnitude = 48._pr*r2inv*r6inv*(r6inv-0.5_pr)
    force_contribution = force_magnitude*particle_separation

    if (measure .and. save_observables) then
        E_potential = E_potential + 4._pr*r6inv*(r6inv-1._pr) - potential_cutoff
        pressure_virial = pressure_virial + particle_distance_squared*force_magnitude
    end if

end subroutine Lennard_Jones

subroutine reaction_field(particle_distance_squared, particle_separation, force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_squared, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff

    real(pr) :: particle_distance, particle_distance_cubed
    real(pr) :: radius_cutoff, radius_cutoff_cubed
    real(pr) :: q_i, q_j, prefac, V_coulomb, V_rf
    real(pr) :: F_coulomb, F_rf, force_magnitude
    real(pr) :: dielectric_factor

    particle_distance = sqrt(particle_distance_squared)
    particle_distance_cubed = particle_distance * particle_distance_squared
    radius_cutoff = sqrt(radius_cutoff_squared)
    radius_cutoff_cubed = radius_cutoff*radius_cutoff_squared

    ! Charges (unit charges)
    q_i = 1._pr
    q_j = 1._pr

    ! Prefactor: 1 / (4*pi*epsilon_0) in atomic units = 1
    prefac = q_i * q_j

    ! Dielectric correction factor (ε → ∞ gives 1)
    dielectric_factor = 1.0_pr

    ! Coulomb energy and force
    V_coulomb = prefac / particle_distance
    F_coulomb = prefac / particle_distance_cubed

    ! Reaction field energy and force correction
    V_rf = prefac * dielectric_factor * (particle_distance_squared - 3._pr*radius_cutoff_squared) / (2._pr * radius_cutoff_cubed)
    F_rf = prefac * dielectric_factor * particle_distance / radius_cutoff_cubed

    ! Total force magnitude
    force_magnitude = F_coulomb + F_rf
    force_contribution = force_magnitude * particle_separation

    if (measure .and. save_observables) then
        E_potential = E_potential + V_coulomb + V_rf - potential_cutoff
        pressure_virial = pressure_virial + particle_distance_squared * force_magnitude
    end if

end subroutine reaction_field

subroutine get_E_potential_contribution(positions, random_particle_id, dE)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_squared
    integer                 :: random_particle_id, i

    dE = 0._pr
    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_squared)
            if (particle_distance_squared <= radius_cutoff_squared) then
                dE = dE + potential_function(particle_distance_squared) ! The term "- potential_cutoff" is irrelevant to the change in potential energy
            end if
        end if
    end do

end subroutine get_E_potential_contribution

END MODULE potentialsMod