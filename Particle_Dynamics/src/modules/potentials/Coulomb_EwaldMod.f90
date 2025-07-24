MODULE Coulomb_EwaldMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod, only: get_distance_squared
    implicit none

    public  Coulomb_realSpace, Coulomb_reciprocalSpace  ! Potential contributions (MC)
    public  Coulomb_Ewald_realSpace, Coulomb_Ewald_reciprocalSpace  ! Full potential
    public  update_Ewald, get_all_expFactors    ! Helper subroutines

    private Coulomb_Ewald_reciprocalSpace_old, get_octant_expFactors    ! Unused/not implemented

CONTAINS

function Coulomb_realSpace(index1, index2, particle_distance_sqr) result(E_potential)
    real(pr), intent(in)   :: particle_distance_sqr
    integer, intent(in)    :: index1, index2
    real(pr)               :: E_potential
    real(pr)               :: particle_distance

    particle_distance = sqrt(particle_distance_sqr)

    E_potential =  ERFC(particle_distance/sigma) * charges(index1) * charges(index2) / particle_distance ! Screened Coulomb term

end function Coulomb_realSpace

function Coulomb_reciprocalSpace(index, proposed_position) result(E_potential)  !! Can be made more efficient
    integer, intent(in)     :: index
    real(pr), intent(in)    :: proposed_position(3)
    real(pr)                :: E_potential
    integer                 :: i
    real(pr)                :: kfactor, kvector(3)
    complex(pr)             :: eikr_old, eikr_new
    complex(pr)             :: kCharge_variation

    E_potential = 0.0_pr

    !$omp parallel private(kvector, kfactor, kCharge_variation, eikr_old, eikr_new) &
    !$omp shared(positions, num_kvec, periodicity, k_vectors, index, charges, reciprocal_charges, proposed_position)&
    !$omp default(none) &
    !$omp reduction(+: E_potential)

    !$omp do schedule(dynamic)
    do i = 1, num_kvec
        kvector = k_vectors(i)%kvector
        kfactor = k_vectors(i)%kfactor

        eikr_old = cmplx(cos(sum(kvector*positions(:,index))), sin(sum(kvector*positions(:,index))), pr)
        eikr_new = cmplx(cos(sum(kvector*proposed_position(:))), sin(sum(kvector*proposed_position(:))), pr)

        kCharge_variation = charges(index) * (eikr_new - eikr_old)

        E_potential  = E_potential + 0.5_pr * kfactor * ( real(kCharge_variation*conjg(kCharge_variation)) + &
            2._pr*real(reciprocal_charges(i)*conjg(kCharge_variation)) )
    end do
    !$omp end do

    !$omp end parallel

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

    term1 = ERFC(particle_distance/sigma) * charges(index1) * charges(index2) / particle_distance ! Screened Coulomb term
    term2 = Ewald_realFactor*exp(-particle_distance_sqr/sigma_sqr) * charges(index1) * charges(index2)

    force_magnitude = term1 + term2
    force_contribution = force_magnitude * particle_separation / particle_distance_sqr

    if (measure .and. save_observables) then
        E_potential = E_potential + term1
        if (debugg .and. (.not. transitory)) E_potential_real = E_potential_real + term1
    end if

end subroutine Coulomb_Ewald_realSpace

subroutine Coulomb_Ewald_reciprocalSpace(force_contribution, E_potential) ! Coulomb potential contribution from other cells in the lattice
    real(pr), intent(out)   :: E_potential, force_contribution(3,num_atoms)
    integer                 :: kx, ky, kz, i, j
    real(pr)                :: factor, kvector(3), kfactor
    complex(pr)             :: eikx(-kgrid(1):kgrid(1), num_atoms)
    complex(pr)             :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr)             :: eikz(-kgrid(3):kgrid(3), num_atoms)
    complex(pr)             :: kCharge, exp_k

    call get_all_expFactors(eikx, eiky, eikz)

    E_potential = 0.0_pr
    force_contribution = 0._pr

    !$omp parallel private(kx, ky, kz, i, j, kCharge, kfactor, exp_k, kvector) &
    !$omp shared(positions, num_atoms, charges, k_vectors, num_kvec, eikx, eiky, eikz, kgrid, k_periodicity) &
    !$omp shared(measure, save_observables) &
    !$omp default(none) &
    !$omp reduction(+: force_contribution, E_potential)

    !$omp do schedule(dynamic)
    do j = 1, num_kvec
        kvector = k_vectors(j)%kvector
        kfactor = k_vectors(j)%kfactor
        kx = k_vectors(j)%kx
        ky = k_vectors(j)%ky
        kz = k_vectors(j)%kz

        kCharge = sum(charges(:)*eikx(kx,:)*eiky(ky,:)*eikz(kz,:))
        do i = 1, num_atoms
            exp_k = eikx(kx,i) * eiky(ky,i) * eikz(kz,i)
            force_contribution(:,i) = force_contribution(:,i) - charges(i) * kfactor * kvector * aimag(kCharge*conjg(exp_k))
        end do
        if (measure .and. save_observables) E_potential  = E_potential + 0.5_pr * kfactor * real(kCharge*conjg(kCharge))
    end do
    !$omp end do

    !$omp end parallel
    if (debugg .and. (.not. transitory)) E_potential_reciprocal = E_potential_reciprocal + E_potential

end subroutine Coulomb_Ewald_reciprocalSpace

subroutine get_all_expFactors(eikx, eiky, eikz)
    complex(pr), intent(out)    :: eikx(-kgrid(1):kgrid(1), num_atoms)
    complex(pr), intent(out)    :: eiky(-kgrid(2):kgrid(2), num_atoms)
    complex(pr), intent(out)    :: eikz(-kgrid(3):kgrid(3), num_atoms)
    integer                     :: i

    ! Calculate exponents with kx, ky, kz = {0, 1}
    eikx(0,:) = (1.0_pr, 0.0_pr)
    eiky(0,:) = (1.0_pr, 0.0_pr)
    eikz(0,:) = (1.0_pr, 0.0_pr)

    ! The shift in positions is done to center the box (code already uses a corner-centered box)
    eikx(1,:) = cmplx(cos(k_periodicity(1)*(positions(1,:))) &
                , sin(k_periodicity(1)*(positions(1,:))), pr)
    eiky(1,:) = cmplx(cos(k_periodicity(2)*(positions(2,:))) &
                , sin(k_periodicity(2)*(positions(2,:))), pr)
    eikz(1,:) = cmplx(cos(k_periodicity(3)*(positions(3,:))) &
                , sin(k_periodicity(3)*(positions(3,:))), pr)

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
    eikx(-kgrid(1):-1,:) = conjg(eikx(kgrid(2):1:-1,:))
    eiky(-kgrid(2):-1,:) = conjg(eiky(kgrid(2):1:-1,:))
    eikz(-kgrid(3):-1,:) = conjg(eikz(kgrid(3):1:-1,:))

end subroutine get_all_expFactors

subroutine update_Ewald(index, proposed_position)
    integer, intent(in)     :: index
    real(pr), intent(in)    :: proposed_position(3)
    real(pr)                :: kvector(3)
    integer                 :: i
    complex(pr)             :: eikr_old, eikr_new

    !$omp parallel private(eikr_old, eikr_new, kvector) &
    !$omp shared(positions, proposed_position, periodicity, index, num_kvec, charges, reciprocal_charges, k_vectors) &
    !$omp default(none)

    !$omp do schedule(dynamic)
    do i = 1, num_kvec
        kvector = k_vectors(i)%kvector

        eikr_old = cmplx(cos(sum(kvector*positions(:,index))), sin(sum(kvector*positions(:,index))), pr)
        eikr_new = cmplx(cos(sum(kvector*proposed_position(:))), sin(sum(kvector*proposed_position(:))), pr)

        reciprocal_charges(i) = reciprocal_charges(i) + charges(index) * (eikr_new - eikr_old)
    end do
    !$omp end do

    !$omp end parallel
end subroutine update_Ewald

!##################################################################################################
!     Not used / Not implemented
!##################################################################################################

subroutine get_octant_expFactors(eikx, eiky, eikz)
    complex(pr), intent(out)    :: eikx(0:kgrid(1), num_atoms)
    complex(pr), intent(out)    :: eiky(0:kgrid(2), num_atoms)
    complex(pr), intent(out)    :: eikz(0:kgrid(3), num_atoms)
    integer                     :: i

    ! Calculate exponents with kx, ky, kz = {0, 1}
    eikx(0,:) = (1.0_pr, 0.0_pr)
    eiky(0,:) = (1.0_pr, 0.0_pr)
    eikz(0,:) = (1.0_pr, 0.0_pr)

    eikx(1,:) = cmplx(cos(k_periodicity(1)*positions(1,:)), sin(k_periodicity(1)*positions(1,:)), pr)
    eiky(1,:) = cmplx(cos(k_periodicity(2)*positions(2,:)), sin(k_periodicity(2)*positions(2,:)), pr)
    eikz(1,:) = cmplx(cos(k_periodicity(3)*positions(3,:)), sin(k_periodicity(3)*positions(3,:)), pr)

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

end subroutine get_octant_expFactors

subroutine Coulomb_Ewald_reciprocalSpace_old(particle_separation, force_contribution, E_potential &
    , pressure_virial)
    real(pr), intent(in)        :: particle_separation(3)
    real(pr), intent(out)       :: force_contribution(3)
    real(pr), intent(inout)     :: E_potential, pressure_virial
    real(pr)                    :: force_magnitude, kx, ky, kz, krx, kry, krz, k_factor
    real(pr)                    :: term1(3), term2(3), term3(3), term4(3)
    integer                     :: i

    do i = 1, num_kvec
        kx = k_vectors(i)%kvector(1)
        ky = k_vectors(i)%kvector(2)
        kz = k_vectors(i)%kvector(3)
        krx = kx * particle_separation(1)
        kry = ky * particle_separation(2)
        krz = kz * particle_separation(3)
        k_factor = k_vectors(i)%kfactor

        ! Every inequivalent contribution (four octants)
        term1 = (/  kx,  ky,  kz/) * sin(  krx + kry + krz)
        if (kx > 0._pr) then
            term2 = (/ -kx,  ky,  kz/) * sin( -krx + kry + krz)
        else
            term2 = 0._pr
        end if
        if (ky > 0._pr) then
            term3 = (/  kx, -ky,  kz/) * sin(  krx - kry + krz)
        else
            term3 = 0._pr
        end if
        if (kz > 0._pr) then
            term4 = (/  kx,  ky, -kz/) * sin(  krx + kry - krz)
        else
            term4 = 0._pr
        end if

        force_contribution = force_contribution + k_factor * (term1 + term2 + term3 + term4)
        if (measure .and. save_observables) E_potential = E_potential + k_factor
    end do

    force_contribution = force_contribution * eightPi_over_volume  ! Accounting for symmetry

    if (measure .and. save_observables) then
        E_potential = twoPi_over_volume * E_potential
    end if


end subroutine Coulomb_Ewald_reciprocalSpace_old

END MODULE Coulomb_EwaldMod