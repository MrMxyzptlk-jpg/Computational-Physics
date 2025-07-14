MODULE potentialContributionMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod
    use potentialPointersMod
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

    dE = E_potential_new - E_potential_old

end subroutine get_E_potential_contribution_normal

subroutine get_E_potential_contribution_Ewald(random_particle_id, proposed_position, dE) ! For MC implementation in long range potentials
    integer, intent(in)     :: random_particle_id
    real(pr), intent(in)    :: proposed_position(3)
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_sqr, E_potential_new, E_potential_old
    integer                 :: i

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

    dE = E_potential_new - E_potential_old

    dE = dE + potential_function_reciprocal()

end subroutine get_E_potential_contribution_Ewald

END MODULE potentialContributionMod