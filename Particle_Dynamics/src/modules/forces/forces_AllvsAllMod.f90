MODULE forces_AllvsAllMod
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

subroutine get_forces_allVSall(E_potential, pressure_virial, pair_corr)
    real(pr), intent(out)   :: E_potential, pressure_virial
    real(pr), intent(inout) :: pair_corr(:)
    real(pr)                :: force_contribution(3), particle_distance_sqr
    integer                 :: i, j

    forces = 0._pr
    if (measure .and. save_observables) then; E_potential = 0.0; pressure_virial = 0.0 ; end if

    !$omp parallel private(j, i, force_contribution, particle_distance_sqr) &
    !$omp shared(positions, num_atoms) &
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

end subroutine get_forces_allVSall

!##################################################################################################
!     Not used / Not implemented
!##################################################################################################

!subroutine get_forces_allVSall_serial(positions, forces, E_potential, pressure_virial, pair_corr) ! Worst performing algorithm
!    real(pr), intent(in)    :: positions(:,:)
!    real(pr), intent(out)   :: forces(:,:)
!    real(pr), intent(out)   :: E_potential, pressure_virial
!    real(pr), intent(inout) :: pair_corr(:)
!    integer(int_huge)       :: i, j
!
!    forces = 0._pr
!    if (measure) then; E_potential = 0.0; pressure_virial = 0.0 ; end if
!
!    do i=1,num_atoms-1
!        do j = i+1, num_atoms
!            call get_force_contribution(i, j, force_contribution, E_potential, &
!                pressure_virial, pair_corr)
!            call add_force_contribution(forces(:,i), forces(:,j), force_contribution, particle_distance_sqr)
!        end do
!    end do
!
!end subroutine get_forces_allVSall_serial

END MODULE forces_AllvsAllMod