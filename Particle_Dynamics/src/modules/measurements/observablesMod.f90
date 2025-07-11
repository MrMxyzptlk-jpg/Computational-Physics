MODULE observablesMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    use meanSquareDisplacementMod
    use pairCorrelationFunctionMod
    use structureFactorMod
    implicit none

CONTAINS

subroutine get_observables(velocities, E_kinetic, Pressure, Temperature)
    real(pr), intent(in)    :: velocities(:,:)
    real(pr), intent(inout) :: Pressure  ! Comes in as Pressure_virial
    real(pr), intent(out)   :: E_kinetic, Temperature

    if (integrator=='velocity-Verlet') then
        E_kinetic = 0.5_pr*sum(velocities*velocities)
        Temperature = Temp_factor*E_kinetic
        Pressure = density*Temperature + Pressure*Pressure_factor
    else
        Pressure = Pressure*Pressure_factor
    end if

end subroutine get_observables

subroutine check_measuring(index, i_measure)
    integer, intent(in)         :: index
    integer, intent(inout)      :: i_measure

    if (mod(index,measuring_jump) == 0 ) then
        measure = .true.
        i_measure = i_measure + 1
    else
        measure = .false.
    end if

end subroutine check_measuring

END MODULE observablesMod