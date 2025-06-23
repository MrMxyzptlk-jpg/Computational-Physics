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

END MODULE subroutinesMod