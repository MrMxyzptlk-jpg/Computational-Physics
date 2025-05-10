MODULE funciones
    USE precision
    implicit none
    
    contains
    
    !funciones de la guía 2
    
    function f1(xx,r)
        real(kind=pr),intent(in),dimension(:)::xx
        real(kind=pr),intent(in)::r
        real(kind=pr),dimension(size(xx))::f1
        f1 = r*xx*(1._pr-xx)
    end function
    
    END module
