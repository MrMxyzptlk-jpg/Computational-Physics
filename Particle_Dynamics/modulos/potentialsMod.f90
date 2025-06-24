MODULE potentialsMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    implicit none

CONTAINS

function Lennard_Jones_potential(particle_distance_sqr) ! Lennard-Jones potential
    real(pr), intent(in)   :: particle_distance_sqr
    real(pr)               :: Lennard_Jones_potential
    real(pr)               :: r2inv, r6inv

    r2inv = 1._pr/particle_distance_sqr
    r6inv = r2inv*r2inv*r2inv
    Lennard_Jones_potential =  4.0_pr *r6inv * (r6inv - 1.0_pr)

end function Lennard_Jones_potential

subroutine Lennard_Jones(particle_distance_sqr, particle_separation,  force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff
    real(pr)                   :: r2inv, r6inv, force_magnitude

    r2inv = 1._pr/particle_distance_sqr
    r6inv = r2inv*r2inv*r2inv
    force_magnitude = 48._pr*r2inv*r6inv*(r6inv-0.5_pr)
    force_contribution = force_magnitude*particle_separation

    if (measure .and. save_observables) then
        E_potential = E_potential + 4._pr*r6inv*(r6inv-1._pr) - potential_cutoff
        pressure_virial = pressure_virial + particle_distance_sqr*force_magnitude
    end if

end subroutine Lennard_Jones

subroutine Coulomb_Ewald_realSpace(particle_distance_sqr, particle_separation, force_contribution, E_potential, pressure_virial) ! Coulomb potential contribution from the reference cell in the lattice
    real(pr), intent(in)        :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)       :: force_contribution(3)
    real(pr), intent(inout)     :: E_potential, pressure_virial
    real(pr)                    :: force_magnitude, particle_distance, term1, term2

    particle_distance = sqrt(particle_distance_sqr)

    term1 = ERFC ( sigma * particle_distance ) / particle_distance ! Screened Coulomb term
    term2 = Ewald_realFactor*exp(-sigma_sqr*particle_distance_sqr)

    force_magnitude = term1 + term2
    force_contribution = force_magnitude * particle_separation / particle_distance_sqr

    if (measure .and. save_observables) then
        E_potential = E_potential + term1
        pressure_virial = pressure_virial + Pressure_factor*particle_distance_sqr*force_magnitude
    end if

end subroutine Coulomb_Ewald_realSpace

subroutine Coulomb_Ewald_reciprocalSpace(particle_distance_sqr, particle_separation, force_contribution, E_potential &
    , pressure_virial) ! Coulomb potential contribution from other cells in the lattice
    real(pr), intent(in)        :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)       :: force_contribution(3)
    real(pr), intent(inout)     :: E_potential, pressure_virial
    real(pr)                    :: force_magnitude, term1(3), term2(3), term3(3), term4(3), krx, kry, krz, k_factor
    integer                     :: i

    do i = 1, num_kvec
        krx = kvectors(i)%kvec(1) * particle_separation(1)
        kry = kvectors(i)%kvec(2) * particle_separation(2)
        krz = kvectors(i)%kvec(3) * particle_separation(3)
        k_factor = kvectors(i)%k_factor

        ! Every inequivalent contribution (four octants)
        term1 = (/  krx,  kry,  krz/) * sin(  krx + kry + krz)
        term2 = (/ -krx,  kry,  krz/) * sin( -krx + kry + krz)
        term3 = (/  krx, -kry,  krz/) * sin(  krx - kry + krz)
        term4 = (/  krx,  kry, -krz/) * sin(  krx + kry - krz)

        force_contribution = force_contribution + k_factor * (term1 + term2 + term3 + term4)
        if (measure .and. save_observables) E_potential = E_potential + k_factor
    end do

    force_contribution = 2._pr * force_contribution * Ewald_forceReciprocalFactor  ! Accounting for symmetry

    if (measure .and. save_observables) then
        E_potential = Ewald_potentialReciprocalFactor * E_potential
        pressure_virial = pressure_virial + Pressure_factor*dot_product(force_contribution, particle_separation)  ! Verify!!
    end if


end subroutine Coulomb_Ewald_reciprocalSpace

subroutine Coulomb_Ewald(particle_distance_sqr, particle_separation,  force_contribution, E_potential, pressure_virial &
    , potential_cutoff)
    real(pr), intent(in)       :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)      :: force_contribution(3)
    real(pr), intent(inout)    :: E_potential, pressure_virial
    real(pr), intent(in)       :: potential_cutoff  ! Irrelevant in Ewald summation, but needed for the interface
    real(pr), dimension(3)     :: force_real, force_reciprocal
    real(pr)                   :: potential_real, potential_reciprocal

    force_real = 0._pr
    force_reciprocal = 0._pr
    potential_real = 0._pr
    potential_reciprocal = 0._pr
    call Coulomb_Ewald_realSpace(particle_distance_sqr, particle_separation, force_real, potential_real, pressure_virial)
    call Coulomb_Ewald_reciprocalSpace(particle_distance_sqr, particle_separation, force_reciprocal, potential_reciprocal &
        , pressure_virial)
    force_contribution = force_real + force_reciprocal

end subroutine Coulomb_Ewald

subroutine get_E_potential_contribution(positions, random_particle_id, dE) ! For MC implementation in short-range potentials
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_sqr
    integer                 :: random_particle_id, i

    dE = 0._pr
    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                dE = dE + potential_function(particle_distance_sqr) ! The term "- potential_cutoff" is irrelevant to the change in potential energy
            end if
        end if
    end do

end subroutine get_E_potential_contribution

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

    ! Prefactor: 1 / (4*pi*epsilon_0) in atomic units = 1
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

END MODULE potentialsMod