MODULE funciones
    USE precision
    implicit none
    
    contains
    
    !funciones de la guía 2

    function hamiltonian_formalism(xx, yy, constants)
        real(kind=pr), intent(in)                        :: xx
        real(kind=pr), dimension(:), intent(in)          :: constants
        real(kind=pr), dimension(:,:), intent(in)        :: yy  ! Input: m x 4
        real(kind=pr), dimension(:,:), allocatable       :: hamiltonian_formalism  ! Output: m x 4
    
        ! Vectorized calculations
        allocate(hamiltonian_formalism(size(yy,1),size(yy,2)))
    
        hamiltonian_formalism(:,1) = yy(:,2)
    
        hamiltonian_formalism(:,2) = -yy(:,1)*(1._pr + 2._pr*constants(1)*yy(:,3)*yy(:,3))
    
        hamiltonian_formalism(:,3) = yy(:,4)
    
        hamiltonian_formalism(:,4) = -yy(:,3)*(1._pr + 2._pr*constants(1)*yy(:,1)*yy(:,1))
    
    end function hamiltonian_formalism

END module
