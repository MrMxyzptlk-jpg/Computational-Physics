MODULE forces_EwaldMod
    use precisionMod
    use constantsMod
    use subroutinesMod
    use parsingMod, only : dim_linkCell, integrator
    use omp_lib
    use variablesMod
    use observablesMod
    use potentialsMod
    use propertiesMod
    use forces_contributionMod
    implicit none

CONTAINS

subroutine get_forces_Ewald(E_potential, pressure_virial, pair_corr)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    integer                 :: i, j
    real(pr)                :: box_dipole(3), surface_potential, pressure_reciprocal
    real(pr)                :: potential_reciprocal
    real(pr)                :: force_contribution(3), particle_distance_sqr
    real(pr)                :: force_reciprocal(3,num_atoms)

    force_reciprocal = 0._pr
    potential_reciprocal = 0._pr

    forces = 0._pr
    if (measure .and. save_observables) then; E_potential = 0.0; pressure_virial = 0.0 ; pressure_reciprocal = 0._pr; end if

    !$omp parallel private(j, i, force_contribution, particle_distance_sqr) &
    !$omp shared(positions, num_atoms) &
    !$omp default(none) &
    !$omp reduction(+: forces, E_potential, pressure_virial, pair_corr)

        !$omp do schedule(dynamic)
        do i=1,num_atoms-1
            do j = i+1, num_atoms
                call get_force_contribution(i, j, force_contribution, E_potential, pressure_virial, pair_corr &
                    , particle_distance_sqr)
                call add_force_contribution(forces(:,i), forces(:,j), force_contribution, particle_distance_sqr)
            end do
        end do
        !$omp end do

    !$omp end parallel

    call Coulomb_Ewald_reciprocalSpace(positions, force_reciprocal, potential_reciprocal)
    forces = forces + force_reciprocal

    ! Correction if there's a net charge (vacuum)
    if(.False.) then ! NOT DEBUGGED AND ONLY EQUIVALENT CHARGES
        forces = forces - pi/(0.75_pr*volume) * spread(sum(positions,2)-0.5_pr*real(num_atoms,pr)*periodicity,2,num_atoms)
        if (measure .and. save_observables) then
            box_dipole = sum(positions,2)-0.5_pr*real(num_atoms,pr)*periodicity ! Net box dipole, when all charges are equal. Else the charges q must multiply each position
            surface_potential  = (pi/(1.5_pr*volume)) * sum(box_dipole*box_dipole)   ! Surface term
            E_potential = E_potential + surface_potential
        end if
    end if

    if (measure .and. save_observables) then
        E_potential = E_potential + potential_reciprocal
        pressure_virial = E_potential
    end if

end subroutine get_forces_Ewald

END MODULE forces_EwaldMod