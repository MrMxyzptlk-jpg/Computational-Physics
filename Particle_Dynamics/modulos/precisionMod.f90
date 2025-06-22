module precisionMod
implicit none
INTEGER, PARAMETER :: int_small             = SELECTED_INT_KIND(2)
INTEGER, PARAMETER :: int_medium            = SELECTED_INT_KIND(4)
INTEGER, PARAMETER :: int_large             = SELECTED_INT_KIND(9)
INTEGER, PARAMETER :: int_huge              = SELECTED_INT_KIND(18)
INTEGER, PARAMETER :: real_simple           = SELECTED_REAL_KIND(p = 6)
INTEGER, PARAMETER :: real_double           = SELECTED_REAL_KIND(p = 15)
INTEGER, PARAMETER :: real_large            = SELECTED_REAL_KIND(p = 30)

integer, parameter :: pr=real_double

END MODULE precisionMod
