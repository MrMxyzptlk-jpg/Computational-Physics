module formats
    implicit none
    character (len=13)   :: format_style = "(*(E14.7,3x))"
    character (len=26)   :: format_style2 = "(2(E14.7,3x),*(E24.16,3x))"
END MODULE formats
