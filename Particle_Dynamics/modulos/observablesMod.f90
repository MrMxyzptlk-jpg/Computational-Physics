MODULE observablesMod
    use precisionMod
    use parametersMod
    use subroutinesMod
    implicit none

    private     positions_buffer, msd_counts, msd_count

    real(pr), allocatable   :: positions_buffer(:,:,:)
    integer, allocatable    :: msd_counts(:)
    integer                 :: msd_count

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

subroutine get_observables(velocities, E_kinetic, Pressure, Temperature)
    real(pr), intent(in)    :: velocities(:,:)
    real(pr), intent(inout) :: Pressure  ! Comes in as Pressure_virial
    real(pr), intent(out)   :: E_kinetic, Temperature

    if (integrator=='velocity-Verlet') then
        E_kinetic = 0.5_pr*sum(velocities*velocities)
        Temperature = Temp_factor*E_kinetic
        Pressure = density*Temperature + Pressure*Pressure_factor
    else
        Pressure = Pressure*Pressure_factor
    end if

end subroutine get_observables

subroutine get_reciprocal_vec(Miller_index, reciprocal_vec)
    integer, intent(in)     :: Miller_index(3)
    real(pr), intent(out)   :: reciprocal_vec(3)
    real(pr)                :: reciprocal_basis(3,3), reciprocal_latticeParam
    integer                 :: i

    select case (structure)
        case ("FCC")
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/-0.5_pr,  0.5_pr,  0.5_pr/)  ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/ 0.5_pr, -0.5_pr,  0.5_pr/)  ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/ 0.5_pr,  0.5_pr, -0.5_pr/)  ! reciprocal lattice basis vector 3
        case ("BCC")
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/0.0_pr, 0.5_pr, 0.5_pr/)      ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/0.5_pr, 0.0_pr, 0.5_pr/)      ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/0.5_pr, 0.5_pr, 0.0_pr/)      ! reciprocal lattice basis vector 3
        case ("random")
            print*, "Random structure selected -> Simple cubic basis set will be employed for structure factor calculation"
            reciprocal_latticeParam = 2._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/1._pr, 0._pr, 0._pr/)         ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/0._pr, 1._pr, 0._pr/)         ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/0._pr, 0._pr, 1._pr/)         ! reciprocal lattice basis vector 3
        case default
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/-0.5_pr,  0.5_pr,  0.5_pr/)   ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/ 0.5_pr, -0.5_pr,  0.5_pr/)   ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/ 0.5_pr,  0.5_pr, -0.5_pr/)   ! reciprocal lattice basis vector 3
    end select

    reciprocal_basis = reciprocal_basis * reciprocal_latticeParam

    reciprocal_vec = 0._pr
    do i =1, 3
        reciprocal_vec = reciprocal_vec + real(Miller_index(i),pr)*reciprocal_basis(:,i)
    end do

end subroutine get_reciprocal_vec

subroutine get_structure_factor(positions, structure_factor, reciprocal_vec)
    real(pr), intent(in)        :: positions(:,:), reciprocal_vec(3)
    real(pr), intent(inout)     :: structure_factor
    complex(pr)                 :: summation
    integer                     :: i

    summation = (0._pr,0._pr)

    do i = 1, num_atoms
        summation = summation + exp(CMPLX(0._pr,sum(reciprocal_vec*positions(:,i)), pr))
    end do

    structure_factor = abs(summation)*abs(summation)/(num_atoms*num_atoms)

end subroutine get_structure_factor

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

subroutine check_measuring(index, i_measure)
    integer, intent(in)         :: index
    integer, intent(inout)      :: i_measure

    if (mod(index,measuring_jump) == 0 ) then
        measure = .true.
        i_measure = i_measure + 1
    else
        measure = .false.
    end if

end subroutine check_measuring

END MODULE observablesMod