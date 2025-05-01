MODULE funciones
    USE precision
    use constantes
    use mzranmod
    use mzranmod_threadsafe
    implicit none

    type :: RNGState
        real(pr) :: next
        logical  :: has_next = .false.
    end type RNGState
    contains

    !funciones de la guía 3

real(pr) function parametrization(x)
    real(kind=pr), intent(in)   ::  x(:)

    if (sum(x*x) > 1._pr) then
        parametrization = 0._pr
    else
        parametrization = 2._pr * sqrt(1._pr - sum(x*x))
    end if

end function

! Returns the analytic volume of a unit n-sphere
real(pr) function volume(n)
    integer, intent(in) :: n

    volume = pi**(real(n, pr)*0.5_pr) / gamma(real(n, pr)*0.5_pr + 1.0_pr)

end function volume

real(pr) function rnd_normal(sigma)
    real(pr), intent(in) :: sigma
    real(pr), save       :: rnd_normal_next
    logical, save        :: has_next = .false.
    real(pr)             :: u1, u2, v1, v2, w, y

    if (has_next) then  ! To avoid wasting the other calculated value
        rnd_normal = sigma * rnd_normal_next
        has_next = .false.
        return
    end if

    do
        u1 = rmzran()
        u2 = rmzran()

        v1 = 2.0_pr * u1 - 1.0_pr
        v2 = 2.0_pr * u2 - 1.0_pr

        w = v1**2 + v2**2
        if (w > 0.0_pr .and. w < 1.0_pr) exit
    end do

    y = sqrt(-2.0_pr * log(w) / w)
    rnd_normal = sigma * v1 * y
    rnd_normal_next = v2 * y
    has_next = .true.
end function rnd_normal

function rnd_normal_parallelized(sigma, state) result(r)
    real(pr), intent(in) :: sigma
    type(MZRanState), intent(inout) :: state
    real(pr) :: u1, u2, v1, v2, w, y
    real(pr), save :: next
    logical, save :: has_next = .false.
    real(pr) :: r

    if (has_next) then
        r = sigma * next
        has_next = .false.
        return
    end if

    do
        u1 = rmzran_threadsafe(state)
        u2 = rmzran_threadsafe(state)

        v1 = 2.0_pr * u1 - 1.0_pr
        v2 = 2.0_pr * u2 - 1.0_pr

        w = v1**2 + v2**2
        if (w > 0.0_pr .and. w < 1.0_pr) exit
    end do

    y = sqrt(-2.0_pr * log(w) / w)
    r = sigma * v1 * y
    next = v2 * y
    has_next = .true.
end function rnd_normal_parallelized

function integrate_trap(lower, upper, N_points, f) result(integral)
    real(pr), intent(in)            :: lower(:), upper(:)
    integer(int_huge), intent(in)   :: N_points(:)
    real(pr)                        :: integral
    integer                         :: dim, i
    integer, allocatable            :: idx(:)
    real(pr), allocatable           :: x(:), h(:)
    real(pr)                        :: weight

    interface
        real(pr)function f(x_x)
            use precision
            implicit none
            real(kind=pr), intent(in)   :: x_x(:)
        end function
    end interface


    dim = size(lower)
    allocate(idx(dim), x(dim), h(dim))

    ! Compute the grid spacings
    h = (upper - lower) / real(N_points-1,pr)

    ! Initialize
    idx = 0
    integral = 0.0_pr

    do
        ! Compute the coordinate vector
        do i = 1, dim
            x(i) = lower(i) + idx(i) * h(i)
        end do

        ! Compute the weight for this point
        weight = 1.0_pr
        do i = 1, dim
            if (idx(i) == 0 .or. idx(i) == N_points(i)-1) then
                weight = weight * 0.5_pr
            end if
        end do

        ! Accumulate
        integral = integral + weight * f(x)

        ! Increment the multi-index
        do i = dim, 1, -1
            idx(i) = idx(i) + 1
            if (idx(i) < N_points(i)) exit
            idx(i) = 0
        end do
        if (all(idx == 0)) exit  ! Full sweep done
    end do

    ! Multiply by the hypervolume element
    integral = integral * product(h)

    deallocate(idx, x, h)
end function integrate_trap

END module