MODULE phase_endMod
    use precisionMod
    use parsingMod
    use checkParsingMod
    use omp_lib
    use writing2filesMod
    use initializationsMod
    use integratorsMod
    use propertiesMod
    use observablesMod
    use measurementsMod
    implicit none

CONTAINS

subroutine phase_end(CPU_elapsed_time)
    real(pr)    :: CPU_elapsed_time

    if (do_pair_correlation) then
        call normalize_pair_correlation(pair_corr)
        call write_pair_corr(pair_corr)
    end if

    if (do_mean_sqr_displacement) then
        call normalize_msd(meanSqrDisplacement)
        call write_msd(meanSqrDisplacement)
    end if

    if (save_state) call write_stateXML()

    call write_output(CPU_elapsed_time, energies(:,1:), pressures(1:), temperatures(1:), structure_factor(1:))

end subroutine phase_end

END MODULE phase_endMod