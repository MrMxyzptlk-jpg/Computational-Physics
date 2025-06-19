MODULE potentialsMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    implicit none

    procedure(pot), pointer :: potential => null()
    procedure(pot_func), pointer :: potential_function => null()

    abstract interface ! Intended to allow for the implementation of a different potential later on
        subroutine pot(particle_distance_squared, particle_separation, force_contribution, E_potential, pressure_virial &
            , potential_cutoff)
            use precisionMod
            real(pr), intent(in)               :: particle_distance_squared, particle_separation(3)
            real(pr), intent(out)              :: force_contribution(3)
            real(pr), intent(inout)            :: E_potential, pressure_virial
            real(pr), intent(in)               :: potential_cutoff
        end subroutine pot
        function pot_func(particle_distance_squared)
            use precisionMod
            real(pr), intent(in)    :: particle_distance_squared
            real(pr)                :: pot_func
        end function pot_func
    end interface

CONTAINS

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

    real(pr) :: r, r2, r3, rcut, rcut3, rcut2
    real(pr) :: q_i, q_j, prefac, V_coulomb, V_rf
    real(pr) :: F_coulomb, F_rf, force_magnitude
    real(pr) :: dielectric_factor

    r2 = particle_distance_squared
    r = sqrt(r2)
    r3 = r2 * r
    rcut = sqrt(radius_cutoff_squared)
    rcut3 = rcut**3
    rcut2 = rcut**2

    ! Charges (unit charges)
    q_i = 1._pr
    q_j = 1._pr

    ! Prefactor: 1 / (4*pi*epsilon_0) in atomic units = 1
    prefac = q_i * q_j

    ! Dielectric correction factor (ε → ∞ gives 1)
    dielectric_factor = 1.0_pr

    ! Coulomb energy and force
    V_coulomb = prefac / r
    F_coulomb = prefac / r3

    ! Reaction field energy and force correction
    V_rf = prefac * dielectric_factor * (r2 - 3*rcut2) / (2 * rcut3)
    F_rf = prefac * dielectric_factor * r / rcut3

    ! Total force magnitude
    force_magnitude = F_coulomb + F_rf
    force_contribution = force_magnitude * particle_separation

    if (measure .and. save_observables) then
        E_potential = E_potential + V_coulomb + V_rf - potential_cutoff
        pressure_virial = pressure_virial + r2 * force_magnitude
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