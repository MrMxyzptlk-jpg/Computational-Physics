! Wrapper module with the different potentials, and a subroutine still to be separated.
MODULE potentialsMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    use subroutinesMod
    use Lennard_JonesMod
    use Coulomb_EwaldMod
    use Coulomb_ReactionFieldMod   ! Not fully implemented yet
    implicit none

CONTAINS

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

END MODULE potentialsMod