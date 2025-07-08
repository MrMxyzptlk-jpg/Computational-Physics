MODULE potentialsMod
    use precisionMod
    use constantsMod
    use parametersMod
    use propertiesMod
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

subroutine Coulomb_Ewald_realSpace(particle_distance_sqr, particle_separation, force_contribution, E_potential, pressure_virial&
    , potential_cutoff) ! Coulomb potential contribution from the reference cell in the lattice
    real(pr), intent(in)        :: particle_distance_sqr, particle_separation(3)
    real(pr), intent(out)       :: force_contribution(3)
    real(pr), intent(inout)     :: E_potential, pressure_virial
    real(pr), intent(in)        :: potential_cutoff  ! Irrelevant in Ewald summation, but needed for the interface
    real(pr)                    :: force_magnitude, particle_distance, term1, term2

    particle_distance = sqrt(particle_distance_sqr)

    term1 = ERFC(particle_distance/sigma) / particle_distance ! Screened Coulomb term
    term2 = Ewald_realFactor*exp(-particle_distance_sqr/sigma_sqr)

    force_magnitude = term1 + term2
    force_contribution = force_magnitude * particle_separation / particle_distance_sqr

    if (measure .and. save_observables) then
        E_potential = E_potential + term1
!        pressure_virial = pressure_virial + Pressure_factor*particle_distance_sqr*force_magnitude
    end if

end subroutine Coulomb_Ewald_realSpace

subroutine Coulomb_Ewald_reciprocalSpace(positions, force_contribution, E_potential) ! Coulomb potential contribution from other cells in the lattice
    real(pr), intent(out)   :: E_potential, force_contribution(3,num_atoms)
    real(pr), intent(in)    :: positions(3,num_atoms)
    integer                 :: kx, ky, kz, k_sqr, i
    real(pr)                :: factor, kvec_real(3), E_contribution
    complex(pr)             :: eikx(        0:kgrid(1), num_atoms)
    complex(pr)             :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr)             :: eikz(-kgrid(3):kgrid(3), num_atoms)
    complex(pr)             :: reciprocal_charge, exp_k
    logical                 :: good_kvec

    ! Calculate exponents with kx, ky, kz = {0, 1}
    eikx(0,:) = (1.0_pr, 0.0_pr)
    eiky(0,:) = (1.0_pr, 0.0_pr)
    eikz(0,:) = (1.0_pr, 0.0_pr)

    eikx(1,:) = cmplx(cos(twoPi*positions(1,:)), sin(twoPi*positions(1,:)), pr)
    eiky(1,:) = cmplx(cos(twoPi*positions(2,:)), sin(twoPi*positions(2,:)), pr)
    eikz(1,:) = cmplx(cos(twoPi*positions(3,:)), sin(twoPi*positions(3,:)), pr)

    ! Use recursion to avoid the calculation of exponential by using multiplication instead
    do i = 2, kgrid(1)
        eikx(i,:) = eikx(i-1,:)*eikx(1,:)
    end do
    do i = 2, kgrid(2)
        eiky(i,:) = eiky(i-1,:)*eiky(1,:)
    end do
    do i = 2, kgrid(3)
        eikz(i,:) = eikz(i-1,:)*eikz(1,:)
    end do

    ! Use conjugation symmetry
    eiky(-kgrid(2):-1,:) = conjg(eiky(kgrid(2):1:-1,:))
    eikz(-kgrid(3):-1,:) = conjg(eikz(kgrid(3):1:-1,:))

    E_potential = 0.0_pr
    force_contribution = 0._pr

    do kx = 0, kgrid(1)
        if (kx == 0) then
            factor = 1.0_pr ! No reflection with respect to y-z plane
        else
            factor = 2.0_pr ! Reflection symmetry with respect to y-z plane
        end if
        ! The same symmetries can be used but might be less efficient (see commented subroutine at the end of this module)
        do ky = -kgrid(2), kgrid(2)
            do kz = -kgrid(3), kgrid(3)
                call check_kvec(kx, ky, kz, k_sqr, good_kvec)
                if (good_kvec) then
                    reciprocal_charge = sum(charges(:)*eikx(kx,:)*eiky(ky,:)*eikz(kz,:))
                    do i = 1, num_atoms
                        exp_k = eikx(kx,i) * eiky(ky,i) * eikz(kz,i)
                        kvec_real = (/kx, ky, kz/) * k_periodicity  ! Real k-space vector
                        force_contribution(:,i) = force_contribution(:,i) &
                            - charges(i) * factor * kfac(k_sqr) * kvec_real * aimag(conjg(reciprocal_charge*exp_k))
                    end do
                    if (measure .and. save_observables) then
                        E_contribution = factor * kfac(k_sqr) * real(reciprocal_charge*conjg(reciprocal_charge))
                        E_potential  = E_potential + E_contribution
                    end if
                end if
            end do
        end do
    end do

    if (measure .and. save_observables) E_potential = E_potential - Ewald_selfTerm

end subroutine Coulomb_Ewald_reciprocalSpace

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

!##################################################################################################
!     Not used / Not implemented
!##################################################################################################

!subroutine Coulomb_Ewald_reciprocalSpace_old(particle_separation, force_contribution, E_potential &
!    , pressure_virial)
!    real(pr), intent(in)        :: particle_separation(3)
!    real(pr), intent(out)       :: force_contribution(3)
!    real(pr), intent(inout)     :: E_potential, pressure_virial
!    real(pr)                    :: force_magnitude, kx, ky, kz, krx, kry, krz, k_factor
!    real(pr)                    :: term1(3), term2(3), term3(3), term4(3)
!    integer                     :: i
!
!    do i = 1, num_kvec
!        kx = kvectors(i)%kvec(1)
!        ky = kvectors(i)%kvec(2)
!        kz = kvectors(i)%kvec(3)
!        krx = kx * particle_separation(1)
!        kry = ky * particle_separation(2)
!        krz = kz * particle_separation(3)
!        k_factor = kvectors(i)%k_factor
!
!        ! Every inequivalent contribution (four octants)
!        term1 = (/  kx,  ky,  kz/) * sin(  krx + kry + krz)
!        if (kx > 0._pr) then
!            term2 = (/ -kx,  ky,  kz/) * sin( -krx + kry + krz)
!        else
!            term2 = 0._pr
!        end if
!        if (ky > 0._pr) then
!            term3 = (/  kx, -ky,  kz/) * sin(  krx - kry + krz)
!        else
!            term3 = 0._pr
!        end if
!        if (kz > 0._pr) then
!            term4 = (/  kx,  ky, -kz/) * sin(  krx + kry - krz)
!        else
!            term4 = 0._pr
!        end if
!
!        force_contribution = force_contribution + k_factor * (term1 + term2 + term3 + term4)
!        if (measure .and. save_observables) E_potential = E_potential + k_factor
!    end do
!
!    force_contribution = force_contribution * eightPi_over_volume  ! Accounting for symmetry
!
!    if (measure .and. save_observables) then
!        E_potential = twoPi_over_volume * E_potential
!    end if
!
!
!end subroutine Coulomb_Ewald_reciprocalSpace_old

END MODULE potentialsMod