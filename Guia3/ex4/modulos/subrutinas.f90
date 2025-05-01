MODULE subrutinas
    use precision
    use omp_lib
    use constantes
    use funciones
    use mzranmod
    use mzranmod_threadsafe
    implicit none

    contains

subroutine MC_volume_gaussian_with_uncertainty_parallel(dim, N, f, estimate, stddev)
    integer, intent(in)          :: dim, N
    real(pr), intent(out)        :: estimate, stddev
    interface
        real(pr) function f(x_x)
            use precision
            real(pr), intent(in) :: x_x(:)
        end function
    end interface

    real(pr)                     :: sigma, gaussian_normalization, gaussian_factor
    real(pr)                     :: rnd_vec(dim - 1)
    real(pr)                     :: rsquare, pdf, fx_over_pdf
    real(pr)                     :: sum_w, sum_w2
    integer                      :: i, j, tid
    type(MZRanState), allocatable :: states(:)
    integer                      :: nthreads

    sigma = 1.0_pr / sqrt(real(dim - 1, pr))
    gaussian_normalization = (1.0_pr / ((2.0_pr * pi * sigma**2)**((dim - 1) / 2.0_pr)))
    gaussian_factor = -1._pr / (2.0_pr * sigma**2)
    sum_w = 0.0_pr
    sum_w2 = 0.0_pr

    ! Determine number of threads and allocate state array
    !$omp parallel
    !$omp single
    nthreads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    allocate(states(nthreads))

    ! Initialize RNG states with unique seeds
    do i = 1, nthreads
        call init_mzran_threadsafe(states(i), 123 + i, 456 + i, 789 + i, 987654321 + i)
    end do

    !$omp parallel default(shared) private(j, i, rsquare, pdf, fx_over_pdf, rnd_vec, tid) &
    !$omp reduction(+:sum_w, sum_w2)
    tid = omp_get_thread_num() + 1

    !$omp do
    do j = 1, N
        rnd_vec = (/(rnd_normal_parallelized(sigma, states(tid)), i=1, dim-1)/)
        rsquare = sum(rnd_vec*rnd_vec)

        if (rsquare <= 1.0_pr) then
            pdf = gaussian_normalization * exp(rsquare * gaussian_factor)
            fx_over_pdf = f(rnd_vec) / pdf
            sum_w  = sum_w  + fx_over_pdf
            sum_w2 = sum_w2 + fx_over_pdf*fx_over_pdf
        end if
    end do
    !$omp end do
    !$omp end parallel

    estimate = sum_w / real(N, pr)
    stddev   = sqrt( (sum_w2 / real(N, pr) - estimate**2) / real(N, pr) )

end subroutine MC_volume_gaussian_with_uncertainty_parallel

subroutine MC_volume_gaussian_with_uncertainty_serial(dim, N, f, estimate, stddev)
    integer, intent(in)          :: dim, N
    real(pr), intent(out)        :: estimate, stddev
    interface
        real(pr) function f(x_x)
            use precision
            real(pr), intent(in) :: x_x(:)
        end function
    end interface

    real(pr)                     :: sigma, gaussian_normalization, gaussian_factor
    real(pr)                     :: rnd_vec(dim - 1)
    real(pr)                     :: rsquare, pdf, fx_over_pdf
    real(pr)                     :: sum_w, sum_w2
    integer                      :: i, j

    sigma = 1.0_pr / sqrt(real(dim - 1, pr))
    gaussian_normalization = (1.0_pr / ((2.0_pr * pi * sigma**2)**((dim - 1) / 2.0_pr)))
    gaussian_factor = -1._pr / (2.0_pr * sigma**2)
    sum_w = 0.0_pr
    sum_w2 = 0.0_pr

    call mzran_init()

    do j = 1, N
        rnd_vec = (/(rnd_normal(sigma), i=1, dim-1)/)
        rsquare = sum(rnd_vec*rnd_vec)

        if (rsquare <= 1.0_pr) then
            pdf = gaussian_normalization * exp(rsquare * gaussian_factor)
            fx_over_pdf = f(rnd_vec) / pdf
            sum_w  = sum_w  + fx_over_pdf
            sum_w2 = sum_w2 + fx_over_pdf*fx_over_pdf
        end if
    end do

    estimate = sum_w / real(N, pr)
    stddev   = sqrt( (sum_w2 / real(N, pr) - estimate**2) / real(N, pr) )

end subroutine MC_volume_gaussian_with_uncertainty_serial

END module subrutinas