MODULE subroutinesMod
    use precisionMod
    use constantsMod
    use parametersMod
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

subroutine update_random_step(N_accepted)
    integer, intent(inout)  :: N_accepted

    if (real(N_accepted,pr)/real(MC_adjust_step,pr) > 0.5_pr) then
        MC_delta = MC_delta*1.05_pr
    else
        MC_delta = MC_delta*0.95_pr
    endif
    N_accepted = 0

end subroutine update_random_step

subroutine update_velocities_velVer(velocities, forces, previous_forces)
    real(pr), dimension(:,:), intent(inout)    :: velocities, forces, previous_forces

    velocities = velocities + dt * 0.5_pr*(previous_forces + forces)

end subroutine update_velocities_velVer

subroutine get_stats(measurements, variance, stddev, average, error)
    real(pr), intent(in)                :: measurements(:)
    real(pr), optional, intent(out)     :: variance, stddev, average, error
    real(pr)                            :: avg, var

    if (present(average) .or. present(variance) .or. present(stddev)) then
        avg = sum(measurements)/real(size(measurements),pr)
    end if

    if (present(variance) .or. present(stddev)) then
        var = sum((measurements - avg)**2) / real(size(measurements) - 1, pr)
    end if

    if (present(average))  average  = avg
    if (present(variance)) variance = var
    if (present(stddev))   stddev   = sqrt(var)
    if (present(error))   error   = sqrt(var/real(size(measurements),pr))

end subroutine get_stats

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

subroutine gasdev_v(harvest)
    real(pr), dimension(:), intent(out) :: harvest
    real(pr), allocatable, save         :: g(:)
    logical, save                       :: gaus_stored = .false.
    integer, save                       :: last_allocated = 0
    real(pr)                            :: rsq, v1, v2
    integer                             :: i, n

    n = size(harvest)

    if (n /= last_allocated) then
        if (last_allocated > 0) deallocate(g)
        allocate(g(n))
        last_allocated = n
        gaus_stored = .false.
    end if

    if (gaus_stored) then
        harvest = g
        gaus_stored = .false.
    else
        i = 1
        do while (i <= n)
            call random_number(v1)
            call random_number(v2)
            v1 = 2.0_pr * v1 - 1.0_pr
            v2 = 2.0_pr * v2 - 1.0_pr
            rsq = v1*v1 + v2*v2

            if (rsq > 0.0_pr .and. rsq < 1.0_pr) then
                rsq = sqrt(-2.0_pr * log(rsq) / rsq)
                harvest(i) = v1 * rsq
                if (i < n) then
                    g(i+1) = v2 * rsq
                end if
                i = i + 2
            end if
        end do
        gaus_stored = .true.
    end if

end subroutine gasdev_v

END MODULE subroutinesMod