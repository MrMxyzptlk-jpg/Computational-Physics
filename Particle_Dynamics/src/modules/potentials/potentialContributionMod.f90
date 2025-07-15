MODULE potentialContributionMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod
    use potentialPointersMod
    use Coulomb_EwaldMod
    implicit none

CONTAINS

subroutine get_E_potential_contribution_normal(random_particle_id, proposed_position,  dE) ! For MC implementation in short-range potentials
    integer, intent(in)     :: random_particle_id
    real(pr), intent(in)    :: proposed_position(3)
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_sqr, E_potential_new, E_potential_old
    integer                 :: i

    E_potential_old = 0._pr
    E_potential_new = 0._pr

    !$omp parallel private(particle_distance_sqr, i) &
    !$omp shared(random_particle_id, num_atoms, positions, radius_cutoff_squared, proposed_position) &
    !$omp reduction(+: E_potential_old, E_potential_new)

    !$omp do schedule(dynamic)
    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                E_potential_old = E_potential_old + potential_function(random_particle_id, i, particle_distance_sqr)
            end if
            call get_distance_squared(proposed_position, positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                E_potential_new = E_potential_new + potential_function(random_particle_id, i, particle_distance_sqr)
            end if
        end if
    end do
    !$omp end do

    !$omp end parallel

    dE = E_potential_new - E_potential_old

end subroutine get_E_potential_contribution_normal

subroutine get_E_potential_contribution_Ewald(random_particle_id, proposed_position, dE) ! For MC implementation in long range potentials
    integer, intent(in)     :: random_particle_id
    real(pr), intent(in)    :: proposed_position(3)
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_sqr, E_potential_new, E_potential_old
    integer                 :: i


    E_potential_old = 0._pr
    E_potential_new = 0._pr

    !$omp parallel private(particle_distance_sqr, i) &
    !$omp shared(random_particle_id, num_atoms, positions, radius_cutoff_squared, proposed_position) &
    !$omp reduction(+: E_potential_old, E_potential_new)

    !$omp do schedule(dynamic)
    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                E_potential_old = E_potential_old + potential_function(random_particle_id, i, particle_distance_sqr)
            end if
            call get_distance_squared(proposed_position, positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                E_potential_new = E_potential_new + potential_function(random_particle_id, i, particle_distance_sqr)
            end if
        end if
    end do
    !$omp end do

    !$omp end parallel

    dE = E_potential_new - E_potential_old

    dE = dE + potential_function_reciprocal(random_particle_id, proposed_position)

end subroutine get_E_potential_contribution_Ewald

subroutine update_potential_contribution_Ewald(index, proposed_position, E_potential, dE)
    real(pr), intent(inout) :: E_potential
    real(pr), intent(in)    :: proposed_position(3), dE
    integer, intent(in)     :: index

    if (measure .and. save_observables) E_potential = E_potential + dE
    MC_accepted = MC_accepted + 1
    positions(:,index) = proposed_position

    call update_Ewald(index, proposed_position)

end subroutine update_potential_contribution_Ewald

subroutine update_potential_contribution_normal(index, proposed_position, E_potential, dE)
    real(pr), intent(inout) :: E_potential
    real(pr), intent(in)    :: proposed_position(3), dE
    integer, intent(in)     :: index

    if (measure .and. save_observables) E_potential = E_potential + dE
    MC_accepted = MC_accepted + 1
    positions(:,index) = proposed_position

end subroutine update_potential_contribution_normal

END MODULE potentialContributionMod