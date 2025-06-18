MODULE integrators
    use precision
    use funciones
    use subrutinas
    use constantes
    use mzranmod
    use parsing
    use linkedLists
    use omp_lib
    use writing2files
    use initializations
    implicit none
contains

subroutine MC_step(i_measure)
    integer(int_large)      :: i_measure

    if (do_linkCell) call create_links(positions)
    call get_forces(positions, forces, Energies(1,i_measure), pressures(i_measure), pair_corr)
    call update_positions_random(positions, Energies(1,i_measure), MC_accepted)

end subroutine MC_step

subroutine velVerlet_step(i_measure)
    integer(int_large)                  :: i_measure

    call update_positions_velVer(positions, velocities, forces)
    previous_forces = forces
    if (do_linkCell) call create_links(positions)
    call get_forces(positions, forces, Energies(1,i_measure), pressures(i_measure), pair_corr)  ! If measure = .False. the observables are ignored

    call update_velocities_velVer(velocities, forces, previous_forces)

end subroutine velVerlet_step

subroutine get_measurements(i_measure)
    integer(int_large)      :: i_measure

    call get_observables(velocities, Energies(2,i_measure), pressures(i_measure), temperatures(i_measure))
    if ((.not.transitory) .and. do_mean_sqr_displacement) call update_msd(positions, meanSqrDisplacement)
    if (do_structure_factor) then
        call get_structure_factor(positions, structure_factor(i_measure), reciprocal_vec)
        call write_tasks(real(i_measure*measuring_jump,pr)*dt, positions, velocities, energies(:,i_measure) &
            , pressures(i_measure), temperatures(i_measure), structure_factor(i_measure))
    else
        call write_tasks(real(i_measure*measuring_jump,pr)*dt, positions, velocities, energies(:,i_measure) &
            , pressures(i_measure), temperatures(i_measure), structure_factor(1))
    end if

end subroutine get_measurements

END MODULE