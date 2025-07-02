program test_epsilon
    implicit none
    real :: x = 3.143
    real(8) :: y = 2.33, epsilon = 0.0
    print *, EPSILON(x)
    print *, EPSILON(y)
end program test_epsilon