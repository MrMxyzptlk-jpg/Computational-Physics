MODULE randomMod
    use precisionMod
    use mzranMod
    use mzran_threadsafeMod
    implicit none

CONTAINS

subroutine gasdev_v(harvest)
    use precisionMod  ! Assumes `pr` is defined here
    implicit none

    real(pr), intent(out)           :: harvest(:)
    real(pr), save                  :: g_saved = 0.0_pr
    logical, save                   :: gaus_stored = .false.
    integer                         :: i, n
    real(pr)                        :: rsq, v1, v2, fac

    n = size(harvest)
    i = 1

    ! Use stored value from last call, if available
    if (gaus_stored) then
        harvest(1) = g_saved
        i = 2
        gaus_stored = .false.
    end if

    ! Generate random pairs using Box-Muller
    do while (i <= n)
        ! Generate two uniform numbers in [-1, 1]
        v1 = 2.0_pr * rmzran() - 1.0_pr
        v2 = 2.0_pr * rmzran() - 1.0_pr
        rsq = v1*v1 + v2*v2

        if (rsq > 0.0_pr .and. rsq < 1.0_pr) then
            fac = sqrt(-2.0_pr * log(rsq) / rsq)
            if (i == n) then
                ! Only one slot left — store one, return the other
                harvest(i) = v1 * fac
                g_saved = v2 * fac
                gaus_stored = .true.
                i = i + 1
            else
                ! Fill two at once
                harvest(i)     = v1 * fac
                harvest(i + 1) = v2 * fac
                i = i + 2
            end if
        end if
    end do

end subroutine gasdev_v

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

END MODULE randomMod