module mzranmod_threadsafe
    use precision
    implicit none

    type :: MZRanState
        integer(int_large) :: i = 521288629
        integer(int_large) :: j = 362436069
        integer(int_large) :: k = 16163801
        integer(int_large) :: n = 1131199299
    end type MZRanState

contains

    function rmzran_threadsafe(state) result(r)
        type(MZRanState), intent(inout) :: state
        real(pr) :: r
        integer(int_large) :: mzran

        mzran = state%i - state%k
        if (mzran < 0) mzran = mzran + 2147483579

        state%i = state%j
        state%j = state%k
        state%k = mzran
        state%n = 69069 * state%n + 1013904243
        mzran = mzran + state%n

        r = real(mzran, kind=pr) * 0.2328306E-9_pr + 0.5_pr
    end function rmzran_threadsafe

    subroutine init_mzran_threadsafe(state, is, js, ks, n)
        type(MZRanState), intent(out) :: state
        integer(int_large), intent(in) :: is, js, ks, n
        state%i = 1 + abs(is)
        state%j = 1 + abs(js)
        state%k = 1 + abs(ks)
        state%n = n
    end subroutine init_mzran_threadsafe

end module mzranmod_threadsafe
