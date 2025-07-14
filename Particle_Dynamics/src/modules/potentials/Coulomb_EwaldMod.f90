MODULE Coulomb_EwaldMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod
    implicit none

    public  Coulomb_Ewald_realSpace, Coulomb_Ewald_reciprocalSpace
    private get_exponential_factors

CONTAINS

function Coulomb_realSpace(index1, index2, particle_distance_sqr) result(E_potential)
    real(pr), intent(in)   :: particle_distance_sqr
    integer, intent(in)    :: index1, index2
    real(pr)               :: E_potential
    real(pr)               :: particle_distance

    particle_distance = sqrt(particle_distance_sqr)

    E_potential =  ERFC(particle_distance/sigma) * charges(index1) * charges(index2) / particle_distance ! Screened Coulomb term

end function Coulomb_realSpace

function Coulomb_reciprocalSpace() result(E_potential)  !! Can be made more efficient
    real(pr)            :: E_potential
    integer             :: kx, ky, kz, k_sqr, i
    real(pr)            :: factor, kvec_real(3)
    complex(pr)         :: eikx(        0:kgrid(1), num_atoms)
    complex(pr)         :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr)         :: eikz(-kgrid(3):kgrid(3), num_atoms)
    complex(pr)         :: reciprocal_charge, exp_k
    logical             :: good_kvec

    call get_exponential_factors(eikx, eiky, eikz)

    E_potential = 0.0_pr

    !$omp parallel private(kx, ky, kz, i, good_kvec, reciprocal_charge, k_sqr, factor) &
    !$omp shared(positions, num_atoms, eikx, eiky, eikz, kfac, charges, kgrid) &
    !$omp reduction(+: E_potential)

    !$omp do schedule(dynamic)
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
                    E_potential  = E_potential + factor * kfac(k_sqr) * real(reciprocal_charge*conjg(reciprocal_charge))
                end if
            end do
        end do
    end do
    !$omp end do

    !$omp end parallel

    E_potential = E_potential - Ewald_selfTerm

end function Coulomb_reciprocalSpace

subroutine Coulomb_Ewald_realSpace(index1, index2, particle_distance_sqr, particle_separation, force_contribution, E_potential &
    , pressure_virial, potential_cutoff) ! Coulomb potential contribution from the reference cell in the lattice
    real(pr), intent(in)    :: particle_distance_sqr, particle_separation(3)
    integer, intent(in)     :: index1, index2
    real(pr), intent(out)   :: force_contribution(3)
    real(pr), intent(inout) :: E_potential, pressure_virial
    real(pr), intent(in)    :: potential_cutoff  ! Irrelevant in Ewald summation, but needed for the interface
    real(pr)                :: force_magnitude, particle_distance, term1, term2

    particle_distance = sqrt(particle_distance_sqr)

    term1 = ERFC(particle_distance/sigma) * charges(index1) * charges(index2)/ particle_distance ! Screened Coulomb term
    term2 = Ewald_realFactor*exp(-particle_distance_sqr/sigma_sqr)

    force_magnitude = term1 + term2
    force_contribution = force_magnitude * particle_separation / particle_distance_sqr

    if (measure .and. save_observables) then
        E_potential = E_potential + term1
    end if

end subroutine Coulomb_Ewald_realSpace

subroutine Coulomb_Ewald_reciprocalSpace(positions, force_contribution, E_potential) ! Coulomb potential contribution from other cells in the lattice
    real(pr), intent(out)   :: E_potential, force_contribution(3,num_atoms)
    real(pr), intent(in)    :: positions(3,num_atoms)
    integer                 :: kx, ky, kz, k_sqr, i
    real(pr)                :: factor, kvec_real(3)
    complex(pr)             :: eikx(        0:kgrid(1), num_atoms)
    complex(pr)             :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr)             :: eikz(-kgrid(3):kgrid(3), num_atoms)
    complex(pr)             :: reciprocal_charge, exp_k
    logical                 :: good_kvec

    call get_exponential_factors(eikx, eiky, eikz)

    E_potential = 0.0_pr
    force_contribution = 0._pr

    !$omp parallel private(kx, ky, kz, i, good_kvec, reciprocal_charge, exp_k, k_sqr, factor, kvec_real) &
    !$omp shared(positions, num_atoms, eikx, eiky, eikz, kgrid, measure, save_observables, kfac, charges, k_periodicity) &
    !$omp reduction(+: force_contribution, E_potential)

    !$omp do schedule(dynamic)
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
                        E_potential  = E_potential + factor * kfac(k_sqr) * real(reciprocal_charge*conjg(reciprocal_charge))
                    end if
                end if
            end do
        end do
    end do
    !$omp end do

    !$omp end parallel

    if (measure .and. save_observables) E_potential = E_potential - Ewald_selfTerm

end subroutine Coulomb_Ewald_reciprocalSpace

subroutine get_exponential_factors(eikx, eiky, eikz)
    complex(pr), intent(out)    :: eikx(        0:kgrid(1), num_atoms)
    complex(pr), intent(out)    :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr), intent(out)    :: eikz(-kgrid(3):kgrid(3), num_atoms)
    integer                     :: i

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

end subroutine get_exponential_factors

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

END MODULE Coulomb_EwaldMod