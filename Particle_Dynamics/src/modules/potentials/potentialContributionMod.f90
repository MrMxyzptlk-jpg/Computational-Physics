MODULE potentialContributionMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod
    use potentialPointersMod
    implicit none

CONTAINS

subroutine get_E_potential_contribution_normal(random_particle_id, dE) ! For MC implementation in short-range potentials
    integer, intent(in)     :: random_particle_id
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_sqr
    integer                 :: i

    dE = 0._pr

    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                dE = dE + potential_function(random_particle_id, i, particle_distance_sqr)
            end if
        end if
    end do

end subroutine get_E_potential_contribution_normal

subroutine get_E_potential_contribution_Ewald(random_particle_id, dE) ! For MC implementation in long range potentials
    integer, intent(in)     :: random_particle_id
    real(pr), intent(out)   :: dE
    real(pr)                :: particle_distance_sqr
    integer                 :: i

    dE = 0._pr

    do i = 1, num_atoms
        if (i /= random_particle_id) then
            call get_distance_squared(positions(:,random_particle_id), positions(:,i), particle_distance_sqr)
            if (particle_distance_sqr <= radius_cutoff_squared) then
                dE = dE + potential_function(random_particle_id, i, particle_distance_sqr)
            end if
        end if
    end do

    dE = dE + potential_function_reciprocal()/num_atoms

end subroutine get_E_potential_contribution_Ewald

END MODULE potentialContributionMod