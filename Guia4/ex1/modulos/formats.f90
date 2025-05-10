module formats
    implicit none
    character (len=13)   :: format_style0 = "(*(E14.7,3x))"
    character (len=22)   :: format_style1 = "((I10,3x),*(E14.7,3x))"
    character (len=33)   :: format_style2 = "((E14.7,3x),(I10,3x),*(E14.7,3x))"

END MODULE formats