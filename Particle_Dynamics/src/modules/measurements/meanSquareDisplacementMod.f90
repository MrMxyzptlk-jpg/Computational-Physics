MODULE meanSquareDisplacementMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    implicit none

    private     positions_buffer, msd_counts, msd_count

    real(pr), allocatable   :: positions_buffer(:,:,:)
    integer, allocatable    :: msd_counts(:)
    integer                 :: msd_count

CONTAINS

subroutine initialize_msd(msd)
    real(pr), allocatable, intent(out)  :: msd(:)

    allocate(positions_buffer(3,num_atoms, 0:max_correlation))
    allocate(msd_counts(0:max_correlation))
    allocate(msd(0:max_correlation))

    positions_buffer = 0._pr
    msd_counts = 0
    msd_count = 0
    msd = 0._pr

end subroutine initialize_msd

subroutine update_msd(positions, msd)
    real(pr), intent(in)    :: positions(:,:)
    real(pr), intent(inout) :: msd(0:)
    integer                 :: i, j, ilast, inext
    real(pr)                :: displacements(3,num_atoms)

    msd_count = msd_count + 1
    ilast = mod(msd_count, max_correlation + 1)

    positions_buffer (:,:,ilast)= positions(:,:)

    do j = 0, min(max_correlation, msd_count-1)
        inext = mod(msd_count - j, max_correlation + 1)

        displacements = positions_buffer(:,:,ilast) - positions_buffer(:,:,inext)
        ! Apply MIC
        do i = 1, 3  ! x, y, z
            displacements(i,:) = displacements(i,:) - periodicity(i) * nint(displacements(i,:) / periodicity(i))
        end do

        msd (j) = msd(j) + sum(displacements*displacements)
        msd_counts(j) = msd_counts(j) + 1
    end do

end subroutine update_msd

subroutine normalize_msd(msd)
    real(pr), intent(out)   :: msd(0:)
    integer                 :: j

    do j = 0, max_correlation
        if (msd_counts(j) > 0) then
            msd(j) = msd(j) / (real(msd_counts(j), pr) * real(num_atoms, pr))
        else
            msd(j) = 0._pr
        end if
    end do

end subroutine normalize_msd

END MODULE meanSquareDisplacementMod