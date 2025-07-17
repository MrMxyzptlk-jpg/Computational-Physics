MODULE subroutinesMod
    use precisionMod
    use constantsMod
    use variablesMod
    use propertiesMod
    implicit none

contains

subroutine get_distance_squared(particle1, particle2, distance_squared)
    real(pr), intent(in)    :: particle1(3), particle2(3)
    real(pr)                :: distance_squared
    real(pr)                :: particle_separation(3)

    particle_separation = particle1 - particle2 ! Separation vector
    particle_separation = particle_separation - periodicity*anint(particle_separation/periodicity) ! PBC
    distance_squared = sum(particle_separation*particle_separation)

end subroutine get_distance_squared

subroutine get_stats(measurements, variance, stddev, average, error)
    real(pr), intent(in)                :: measurements(:)
    real(pr), optional, intent(out)     :: variance, stddev, average, error
    real(pr)                            :: avg, var

    if (.not. any((/present(average), present(variance), present(stddev), present(error)/))) then
        print*, "Unspecified stat inquiry"
        return
    end if

    avg = sum(measurements)/real(size(measurements),pr)

    if (present(average))  average  = avg

    if (any((/present(variance), present(stddev), present(error)/))) then
        var = sum((measurements - avg)**2) / real(size(measurements) - 1, pr)
        if (present(variance)) variance = var
        if (present(stddev))   stddev   = sqrt(var)
        if (present(error))   error   = sqrt(var/real(size(measurements),pr))
    end if


end subroutine get_stats

subroutine check_kvec(nx, ny, nz, k_sqr, kvec_flag) ! Verify the k-vector is inside the 3D ellipsoid in reciprocal space
    integer, intent(in)     :: nx, ny, nz
    real(pr), intent(out)   :: k_sqr
    logical, intent(out)    :: kvec_flag
    real(pr)                :: kx, ky, kz

    kvec_flag = .False.

    kx = real(nx*kgrid(1),pr)
    ky = real(ny*kgrid(2),pr)
    kz = real(nz*kgrid(3),pr)
    k_sqr = kx*kx + ky*ky + kz*kz

    if ((k_sqr <= 1._pr) .and. (k_sqr /= 0._pr)) kvec_flag = .True.

end subroutine check_kvec

END MODULE subroutinesMod