MODULE potentialPointersMod
    use precisionMod
    implicit none

    ! Potentials' pointers
    procedure(pot), pointer      :: potential => null()
    procedure(pot_func), pointer :: potential_function => null()
    procedure(pot_func_recip), pointer  :: potential_function_reciprocal => null()
    procedure(pot_contrib), pointer     :: get_E_potential_contribution => null()

    abstract interface ! Intended to allow for the implementation of a different potential later on
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

        function pot_func_recip()
            use precisionMod
            real(pr)                :: pot_func_recip
        end function pot_func_recip

        subroutine pot_contrib(random_particle_id, dE)
            use precisionMod
            integer, intent(in)     :: random_particle_id
            real(pr), intent(out)   :: dE
        end subroutine pot_contrib
    end interface

END MODULE potentialPointersMod