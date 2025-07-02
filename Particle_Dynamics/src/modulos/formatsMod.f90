module formatsMod
    implicit none
    character (len=13)   :: format_style0 = "(*(E14.7,3x))"
    character (len=16)   :: format_XYZ   = "(A5,*(E14.7,3x))"
    character (len=26)   :: format_style2 = "(2(E14.7,3x),*(E24.16,3x))"
    character (len=18)   :: format_observables = "(a,*(a,E14.8,3x,))"
    character (len=21)   :: format_state = "(2(E24.16,3x),E24.16)"
END MODULE formatsMod
