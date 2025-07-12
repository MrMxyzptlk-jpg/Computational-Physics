MODULE potentialContributionMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod
    implicit none

CONTAINS

subroutine get_E_potential_contribution_normal(random_particle_id, dE) ! For MC implementation in short-range potentials
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

    dE = dE * 0.5_pr  ! Avoid double-counting

end subroutine get_E_potential_contribution_normal

subroutine get_E_potential_contribution_Ewald(random_particle_id, dE) ! For MC implementation in short-range potentials
    use Coulomb_EwaldMod
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

    call get_exponential_factors

    dE = dE * 0.5_pr  ! Avoid double-counting

end subroutine get_E_potential_contribution_Ewald

END MODULE potentialContributionMod