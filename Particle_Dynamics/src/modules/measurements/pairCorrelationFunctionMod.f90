MODULE pairCorrelationFunctionMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    implicit none

CONTAINS

subroutine get_pair_correlation(positions, pair_corr)
    real(pr), intent(in)        :: positions(:,:)
    real(pr), intent(inout)     :: pair_corr(:)
    real(pr)                    :: particle_distance_squared
    integer(int_large)          :: i,j

    if (measure) then
        !$omp parallel private(j, i,  particle_distance_squared) &
        !$omp shared(positions, num_atoms) &
        !$omp reduction(+: pair_corr)

            !$omp do schedule(dynamic)
            do i=1,num_atoms-1
                do j = i+1, num_atoms
                    call get_distance_squared(positions(:,i), positions(:,j), particle_distance_squared)
                    call update_pair_correlation(particle_distance_squared, pair_corr)
                end do
            end do
            !$omp end do

        !$omp end parallel
    end if

end subroutine get_pair_correlation

subroutine update_pair_correlation(particle_distance_squared, pair_corr)
    real(pr), intent(in)        :: particle_distance_squared
    real(pr), intent(inout)     :: pair_corr(:)
    integer(int_large)          :: radius_bin

    if (measure) then
        if (particle_distance_squared <= pair_corr_cutoff_sqr) then
            radius_bin = int(sqrt(particle_distance_squared) / dr) + 1
            pair_corr(radius_bin) = pair_corr(radius_bin) + 1._pr
        end if
    end if

end subroutine update_pair_correlation

subroutine normalize_pair_correlation(pair_corr)
    real(pr), intent(inout)     :: pair_corr(:)
    integer                     :: i
    real(pr)                    :: r_lower, r_upper, shell_vol, ideal_pair_corr

    do i = 1, pair_corr_bins
        r_lower = real(i-1,pr) * dr
        r_upper = r_lower + dr
        shell_vol = (4.0_pr*pi / 3.0_pr) * (r_upper**3 - r_lower**3)
        ideal_pair_corr = density * shell_vol * real(num_atoms,pr)
        pair_corr(i) = pair_corr(i) * 2._pr / ideal_pair_corr / real(measuring_steps,pr)
    end do

end subroutine normalize_pair_correlation

END MODULE pairCorrelationFunctionMod