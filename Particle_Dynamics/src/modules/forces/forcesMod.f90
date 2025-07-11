! Wrapper module for different force summations
MODULE forcesMod
    use precisionMod
    use forces_EwaldMod
    use forces_LinkedListsMod
    use forces_AllvsAllMod
    implicit none

    abstract interface

        subroutine force_sub(E_potential, pressure_virial, pair_corr)
            use precisionMod
            real(pr), intent(out)   :: E_potential, pressure_virial
            real(pr), intent(inout) :: pair_corr(:)
        end subroutine force_sub
    end interface

    procedure(force_sub), pointer   :: get_forces          => null()

END MODULE forcesMod