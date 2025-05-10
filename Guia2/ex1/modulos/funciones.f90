MODULE funciones
    USE precision
    USE constantes
    implicit none

    contains

    !funciones de la guía 2

    function analitical_solution(x,t,num_sum_term, initial_Temp) ! Analitical solution to the non-dimensional problem
        real(kind=pr), dimension(:), intent(in)     :: x
        real(kind=pr), dimension(size(x))           :: analitical_solution
        real(kind=pr), intent(in)                   :: t, initial_Temp
        real(kind=pr)                               :: common_factor,K_n
        integer(kind=int_huge), intent(in)          :: num_sum_term
        integer(kind=int_huge)                      :: i

        analitical_solution = 0._pr
        common_factor = initial_Temp*4._pr/pi

        do i = 2*num_sum_term-1, 1, -2 ! Sum in inverse order to avoid precision errors. Skip end points due to the boundary conditions
            K_n = pi*real(i,pr)
            analitical_solution(2:size(x)-1)= analitical_solution(2:size(x)-1) + (1._pr/real(i,pr))&
            *sin(pi*real(i,pr)*x(2:size(x)-1))*exp(-K_n*K_n*t)
        end do

        analitical_solution = analitical_solution*common_factor

    end function analitical_solution

    function analitical_solution_term(x,t,num_sum_term, initial_Temp) ! Analitical solution to the non-dimensional problem
        real(kind=pr), dimension(:), intent(in)     :: x
        real(kind=pr), dimension(size(x))           :: analitical_solution_term
        real(kind=pr), intent(in)                   :: t, initial_Temp
        real(kind=pr)                               :: common_factor,K_n
        integer(kind=int_huge), intent(in)          :: num_sum_term

        analitical_solution_term = 0._pr
        common_factor = initial_Temp*4._pr/pi

        K_n = pi*real(num_sum_term,pr)
        analitical_solution_term = analitical_solution_term + (1._pr/real(num_sum_term,pr))*sin(K_n*x)*exp(-K_n*K_n*t)

        analitical_solution_term = analitical_solution_term*common_factor

    end function analitical_solution_term

    END module