MODULE funciones
    USE precision
    implicit none
    
    contains
    
    !funciones de la guía 2

    function lagrangian_formalism(xx, yy, constants)
        real(kind=pr), intent(in)                        :: xx
        real(kind=pr), dimension(:), intent(in)          :: constants
        real(kind=pr), dimension(:,:), intent(in)        :: yy  ! Input: m x 4
        real(kind=pr), dimension(:,:), allocatable       :: lagrangian_formalism  ! Output: m x 4
        real(kind=pr), dimension(size(yy,1))             :: sin_diff, cos_diff, denominator
        
        ! Vectorized calculations
        allocate(lagrangian_formalism(size(yy,1),size(yy,2)))
        
        sin_diff = sin(yy(:,1) - yy(:,3))
        cos_diff = cos(yy(:,1) - yy(:,3))
        denominator = 1._pr + constants(1) * sin_diff * sin_diff
    
        lagrangian_formalism(:,1) = yy(:,2)
    
        lagrangian_formalism(:,2) = (-(1._pr + constants(1)) * constants(3) * sin(yy(:,1)) - constants(1) * constants(2) &
        *yy(:,4) * yy(:,4) * sin_diff - constants(1) * cos_diff * (yy(:,2) * yy(:,2) * sin_diff - constants(3) * sin(yy(:,3)))) &
        / denominator
    
        lagrangian_formalism(:,3) = yy(:,4)
    
        lagrangian_formalism(:,4) = ((1._pr + constants(1)) * (yy(:,2) * yy(:,2) * sin_diff - constants(3) * sin(yy(:,3))) &
        + cos_diff * ((1._pr + constants(1)) * constants(3) * sin(yy(:,1)) + constants(1) * constants(2) * yy(:,4) * yy(:,4) &
        * sin_diff)) / (constants(2) * denominator)
    
    end function lagrangian_formalism

    function compute_p1(yy, constants)
        real(kind=pr), dimension(:), intent(in)         :: constants
        real(kind=pr), dimension(:,:), intent(in)       :: yy  ! Input: m x 4
        real(kind=pr), dimension(size(yy,1))            :: compute_p1 ! Output: m x 4

        compute_p1 = (1._pr + constants(1))*yy(:,2) + constants(1)*constants(2)*yy(:,4)*cos(yy(:,1) - yy(:,3))

    end function compute_p1

    function compute_p2(yy, constants)
        real(kind=pr), dimension(:), intent(in)         :: constants
        real(kind=pr), dimension(:,:), intent(in)       :: yy  ! Input: m x 4
        real(kind=pr), dimension(size(yy,1))            :: compute_p2 ! Output: m x 4

        compute_p2 = constants(1)*constants(2)*(constants(2)*yy(:,4) + yy(:,2)*cos(yy(:,1) - yy(:,3)))

    end function compute_p2

!############################################################################################################## 
! NOT USED:
!############################################################################################################## 

    function hamiltonian_formalism(xx, yy, constants)
        real(kind=pr), intent(in)                        :: xx
        real(kind=pr), dimension(:), intent(in)          :: constants
        real(kind=pr), dimension(:,:), intent(in)        :: yy  ! Input: m x 4
        real(kind=pr), dimension(:,:), allocatable       :: hamiltonian_formalism  ! Output: m x 4
        real(kind=pr), dimension(size(yy,1))             :: sin_diff, cos_diff, denominator, h1, h2
    
        ! Vectorized calculations
        allocate(hamiltonian_formalism(size(yy,1),size(yy,2)))

        
        sin_diff = sin(yy(:,1) - yy(:,3))
        cos_diff = cos(yy(:,1) - yy(:,3))
        denominator = (1._pr + constants(1) * sin_diff*sin_diff)*constants(2)
        h1 = yy(:,2)*yy(:,4)*sin_diff/denominator
        h2 = (yy(:,2)*yy(:,2)*constants(1)*constants(2)*constants(2) + (1._pr + constants(1))*yy(:,4)*yy(:,4) - 2._pr &
        *constants(1)*constants(2)*yy(:,2)*yy(:,4)*cos_diff ) / (2._pr*denominator*denominator)

        hamiltonian_formalism(:,1) = (constants(2)*yy(:,2)-yy(:,4)*cos_diff)/denominator
    
        hamiltonian_formalism(:,2) = -constants(3)*(1._pr + constants(1))*sin(yy(:,1)) - h1 + h2* sin(2._pr*(yy(:,1) - yy(:,3)))
    
        hamiltonian_formalism(:,3) = (-constants(1)*constants(2)*yy(:,2)*cos_diff + yy(:,4)*(1._pr + constants(1))) &
        /(denominator*constants(1))
    
        hamiltonian_formalism(:,4) = -constants(3)*constants(2)*constants(1)*sin(yy(:,3)) + h1 - h2* sin(2._pr*(yy(:,1) - yy(:,3)))
    
    end function hamiltonian_formalism

    function compute_theta1(yy, constants)
        real(kind=pr), dimension(:), intent(in)         :: constants
        real(kind=pr), dimension(:,:), intent(in)       :: yy  ! Input: m x 4
        real(kind=pr), dimension(size(yy,1))            :: compute_theta1 ! Output: m x 4

        compute_theta1 = (constants(2)*yy(:,2) - yy(:,4)*cos(yy(:,1)-yy(:,3))) / (constants(2)*(1._pr + constants(1)&
        *sin(yy(:,1)-yy(:,3))*sin(yy(:,1)-yy(:,3))))

    end function compute_theta1

    function compute_theta2(yy, constants)
        real(kind=pr), dimension(:), intent(in)         :: constants
        real(kind=pr), dimension(:,:), intent(in)       :: yy  ! Input: m x 4
        real(kind=pr), dimension(size(yy,1))            :: compute_theta2 ! Output: m x 4

        compute_theta2 = (-constants(2)*constants(1)*yy(:,2)*cos(yy(:,1)-yy(:,3)) + yy(:,4)*(1._pr + constants(1))) &
        / (constants(2)*constants(2)*constants(1)*(1._pr + constants(1)*sin(yy(:,1)-yy(:,3))*sin(yy(:,1)-yy(:,3))))

    end function compute_theta2

    END module
