MODULE potentialPointersMod
    use precisionMod
    implicit none

    ! Potentials' pointers
    procedure(pot), pointer             :: potential => null()  ! To get the full system potential energy
    procedure(pot_func), pointer        :: potential_function => null() ! To get short range contribution to the potential difference in MC
    procedure(pot_func_recip), pointer  :: potential_function_reciprocal => null()  ! To get the long range contribution to the potential difference in MC
    procedure(get_potContrib), pointer  :: get_E_potential_contribution  => null()  ! To get the potential energy contribution for an arbitrary potential in MC trial
    procedure(acceptTrial), pointer   :: update_potential_contribution => null()  ! To accept a MC trial

    abstract interface ! Intended to allow for the implementation of a different potentials later on
        subroutine pot(index1, index2, particle_distance_squared, particle_separation, force_contribution, E_potential &
            , pressure_virial, potential_cutoff)
            use precisionMod
            real(pr), intent(in)    :: particle_distance_squared, particle_separation(3)
            real(pr), intent(out)   :: force_contribution(3)
            real(pr), intent(inout) :: E_potential, pressure_virial
            real(pr), intent(in)    :: potential_cutoff
            integer, intent(in)     :: index1, index2
        end subroutine pot

        function pot_func(index1, index2, particle_distance_squared)
            use precisionMod
            real(pr), intent(in)    :: particle_distance_squared
            integer, intent(in)     :: index1, index2
            real(pr)                :: pot_func
        end function pot_func

        function pot_func_recip(random_particle_id, proposed_position)
            use precisionMod
            integer, intent(in)     :: random_particle_id
            real(pr), intent(in)    :: proposed_position(3)
            real(pr)                :: pot_func_recip
        end function pot_func_recip

        subroutine get_potContrib(random_particle_id, proposed_position, dE)
            use precisionMod
            integer, intent(in)     :: random_particle_id
            real(pr), intent(in)    :: proposed_position(3)
            real(pr), intent(out)   :: dE
        end subroutine get_potContrib

        subroutine acceptTrial(index, proposed_position, E_potential, dE)
            use precisionMod
            real(pr), intent(inout) :: E_potential
            real(pr), intent(in)    :: proposed_position(3), dE
            integer, intent(in)     :: index
        end subroutine acceptTrial
    end interface

END MODULE potentialPointersMod