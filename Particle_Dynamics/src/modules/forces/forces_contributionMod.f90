MODULE forces_contributionMod
    use precisionMod
    use constantsMod
    use subroutinesMod
    use parsingMod, only : dim_linkCell, integrator
    use omp_lib
    use variablesMod
    use observablesMod
    use potentialsMod
    use propertiesMod
    implicit none

CONTAINS

subroutine get_force_contribution(index1, index2, force_contribution, E_potential, pressure_virial, pair_corr &
    , particle_distance_sqr)
    integer, intent(in)                     :: index1, index2
    real(pr), dimension(3), intent(out)     :: force_contribution
    real(pr), intent(out)                   :: E_potential, pressure_virial, particle_distance_sqr
    real(pr), intent(inout)                 :: pair_corr(:)
    real(pr), dimension(3)                  :: particle_separation

    particle_separation = positions(:,index1) - positions(:,index2) ! Separation vector
    particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
    particle_distance_sqr = sum(particle_separation*particle_separation)

    if (particle_distance_sqr <= radius_cutoff_squared) then
        call potential(index1, index2, particle_distance_sqr, particle_separation, force_contribution, E_potential &
            , pressure_virial, potential_cutoff)
    endif

    if (do_pair_correlation .and. .not. transitory) call update_pair_correlation(particle_distance_sqr, pair_corr)

end subroutine get_force_contribution

subroutine add_force_contribution(particle1_forces, particle2_forces, force_contribution, particle_distance_sqr)
    real(pr), dimension(3), intent(out)     :: particle1_forces, particle2_forces
    real(pr), dimension(3), intent(in)      :: force_contribution
    real(pr), intent(in)                    :: particle_distance_sqr


    if (particle_distance_sqr <= radius_cutoff_squared) then
        particle1_forces = particle1_forces + force_contribution
        particle2_forces = particle2_forces - force_contribution
    endif

end subroutine add_force_contribution

END MODULE forces_contributionMod