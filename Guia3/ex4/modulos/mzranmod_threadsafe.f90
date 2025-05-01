module mzranmod_threadsafe
    implicit none
    integer, parameter :: K4B = selected_int_kind(9)
    integer, parameter :: DP  = kind(1.0d0)

    type :: MZRanState
        integer(K4B) :: i = 521288629
        integer(K4B) :: j = 362436069
        integer(K4B) :: k = 16163801
        integer(K4B) :: n = 1131199299
    end type MZRanState

contains

    function rmzran_threadsafe(state) result(r)
        type(MZRanState), intent(inout) :: state
        real(DP) :: r
        integer(K4B) :: mzran

        mzran = state%i - state%k
        if (mzran < 0) mzran = mzran + 2147483579

        state%i = state%j
        state%j = state%k
        state%k = mzran
        state%n = 69069 * state%n + 1013904243
        mzran = mzran + state%n

        r = real(mzran, kind=DP) * 0.2328306E-9_DP + 0.5_DP
    end function rmzran_threadsafe

    subroutine init_mzran_threadsafe(state, is, js, ks, n)
        type(MZRanState), intent(out) :: state
        integer(K4B), intent(in) :: is, js, ks, n
        state%i = 1 + abs(is)
        state%j = 1 + abs(js)
        state%k = 1 + abs(ks)
        state%n = n
    end subroutine init_mzran_threadsafe

end module mzranmod_threadsafe
