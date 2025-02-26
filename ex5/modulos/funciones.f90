MODULE funciones
    USE ISO
    USE precision
    implicit none
    
    contains
    
    !funciones de la guía 5
    
        function f1(xx,yy)
            real(kind=pr),intent(in)::xx, yy
            real(kind=pr)::f1
            f1=-xx*yy
        end function
    
        function f2(xx)
            real(kind=pr),intent(in)::xx
            real(kind=pr)::f2
            f2=exp(-xx**2._pr/2._pr)
        end function
    
    END module