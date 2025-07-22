MODULE dimensionsMod
    use precisionMod
    use variablesMod
    implicit none

CONTAINS

elemental real(pr) function redimensionalize(x, quantity)
    real(pr), intent(in)            :: x
    character(len=*), intent(in)    :: quantity

    select case(quantity)
        case ("distance")
            redimensionalize = x * conversion_factors(1)    ! Distance      [Bohr = a₀]
        case ("time")
            redimensionalize = x * conversion_factors(2)    ! Time          [ℏ/Eh = tₐ]
        case ("temperature")
            redimensionalize = x * conversion_factors(3)    ! Temperature   [Eh / kB]
        case ("energy")
            redimensionalize = x * conversion_factors(4)    ! Energy        [hartree = Eh]
        case ("velocity")
            redimensionalize = x * conversion_factors(5)    ! Velocity      [a₀ / tₐ]
        case ("pressure")
            redimensionalize = x * conversion_factors(6)    ! Pressure      [hartree / bohr³]
        case ("volume")
            redimensionalize = x * conversion_factors(7)    ! Volume        [Bohr³]
        case ("density")
            redimensionalize = x * conversion_factors(8)    ! Density       [1/Bohr³]
    end select

end function redimensionalize

elemental real(pr) function adimensionalize(x, quantity)
    real(pr), intent(in)            :: x
    character(len=*), intent(in)    :: quantity

    select case(quantity)
        case ("distance")
            adimensionalize = x / conversion_factors(1)     ! Distance      [Bohr = a₀]
        case ("time")
            adimensionalize = x / conversion_factors(2)     ! Time          [ℏ/Eh = tₐ]
        case ("temperature")
            adimensionalize = x / conversion_factors(3)     ! Temperature   [Eh / kB]
        case ("energy")
            adimensionalize = x / conversion_factors(4)     ! Energy        [hartree = Eh]
        case ("velocity")
            adimensionalize = x / conversion_factors(5)     ! Velocity      [a₀ / tₐ]
        case ("pressure")
            adimensionalize = x / conversion_factors(6)     ! Pressure      [hartree / bohr³]
        case ("volume")
            adimensionalize = x / conversion_factors(7)     ! Volume        [Bohr³]
        case ("density")
            adimensionalize = x / conversion_factors(8)     ! Density       [1/Bohr³]
    end select

end function adimensionalize

END MODULE dimensionsMod