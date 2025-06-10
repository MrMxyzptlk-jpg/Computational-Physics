MODULE writing2files
    use precision
    use subrutinas
    use parsing
    implicit none

    private     size_x, size_y, size_z, symbol
    character(len=11)       :: size_x, size_y, size_z
    character(len=3)        :: symbol

contains

subroutine initialize_XYZ_data()

    write (size_x,'(E11.5)') periodicity(1)
    write (size_y,'(E11.5)') periodicity(2)
    write (size_z,'(E11.5)') periodicity(3)

    size_x = adjustl(trim(size_x))
    size_y = adjustl(trim(size_y))
    size_z = adjustl(trim(size_z))

end subroutine initialize_XYZ_data

subroutine write_XYZfile(positions, velocities, time, unitnum)
    real (pr), intent (in)              :: positions(:,:), velocities(:,:), time
    integer(int_medium), intent (in)    :: unitnum
    integer                             :: i
    character(len=11)                   :: time_tmp

    symbol = 'X'      ! Change to real element if needed

    ! Write time in string format
    write(time_tmp,'(E11.5)') time*conversion_factors(2)
    time_tmp = adjustl(trim(time_tmp))

    ! Line 1: number of particles
    write(unitnum, '(i6)') num_atoms

    ! Line 2: extended XYZ header with box info and time
    write(unitnum,'(A)') 'Lattice="' // &
         size_x // ' 0.0  0.0  0.0 ' // &
         size_y // ' 0.0  0.0  0.0 ' // &
         size_z // '" Properties=species:S:1:pos:R:3:vel:R:3 Time=' // &
         time_tmp


    do i = 1, num_atoms
        write(unitnum, fmt=format_XYZ) symbol, positions(:,i)*conversion_factors(1), velocities(:,i)*conversion_factors(5)
    end do

end subroutine write_XYZfile

subroutine write_output(CPU_elapsed_time, energies, pressures, temperatures)
    real (pr), intent (in)  :: CPU_elapsed_time, energies(:,:), pressures(:), temperatures(:)
    real (pr)               :: energy_avg(2), pressure_avg, temperature_avg
    real (pr)               :: energy_stddev(2), pressure_stddev, temperature_stddev
    integer                 :: unit_info

    open(newunit=unit_info, file="datos/INFO.out", status="replace")
        write(unit_info,'(a20,13x,I4)')     "Number of atoms:       ", num_atoms
        write(unit_info,'(a20,14x,a)')      "Initial structure:     ", structure
        write(unit_info,'(a20,4x,a)')       "Potential:             ", type
        write(unit_info,'(a20,2x,a)')       "Integrator:            ", integrator
        write(unit_info,'(a20,11x,I6)')     "Transitory steps:      ", transitory_steps
        write(unit_info,'(a20,11x,I6)')     "Run steps:             ", MD_steps
        write(unit_info,'(a20,6x,E11.5)')   "Simulated time:        ", real(MD_steps,pr)*dt*conversion_factors(2)

        if (save_observables) then
            call get_stats(energies(1,:), average = energy_avg(1), stddev = energy_stddev(1))
            write(unit_info,format_observables)     "Potential Energy:  ", " Average = ", energy_avg(1)*conversion_factors(4)  &
                , " Standard deviation = ", energy_stddev(1)*conversion_factors(4)
            call get_stats(energies(2,:), average = energy_avg(2), stddev = energy_stddev(2))
            write(unit_info,format_observables)     "Kinetic Energy:    ", " Average = ", energy_avg(2)*conversion_factors(4)  &
                , " Standard deviation = ", energy_stddev(2)*conversion_factors(4)
            call get_stats(pressures, average = pressure_avg, stddev = pressure_stddev)
            write(unit_info,format_observables)     "Pressure:          ", " Average = ", pressure_avg*conversion_factors(6)  &
                , " Standard deviation = ", pressure_stddev*conversion_factors(6)
            call get_stats(temperatures, average = temperature_avg, stddev = temperature_stddev)
            write(unit_info,format_observables)     "Temperature:       ", " Average = ", temperature_avg*conversion_factors(3)  &
                , " Standard deviation = ", temperature_stddev*conversion_factors(3)
        end if

        write(unit_info,'(a)')"-----------------------------------------------------------------------------"
        write(unit_info,'(a20,5x,F11.5, a)')"Elapsed time:          ", CPU_elapsed_time,"s"
    close(unit_info)

end subroutine write_output

subroutine write_pair_corr(pair_corr)
    real (pr), intent (in)  :: pair_corr(:)
    real (pr)               :: bin_center
    integer                 :: unitnum, i

    open(newunit=unitnum, file="datos/pair_correlation.out", status="replace")
        write(unitnum, '(a)') "## r | pair_correlation(r)"
        do i =1, pair_corr_bins
            bin_center = (real(i-1,pr) + 0.5_pr)*dr
            write(unitnum, format_style0) bin_center, pair_corr(i)
        end do
    close(unitnum)

end subroutine write_pair_corr


subroutine write_structure_factor(structure_factor, reciprocal_vec)
    real (pr), intent (in)  :: structure_factor(:), reciprocal_vec(3)
    integer                 :: unitnum, i

    open(newunit=unitnum, file="datos/structure_factor.out", status="replace")
        write(unitnum, '(a,3(E12.5,1x))') "## Reciprocal vector: K = ", reciprocal_vec
        write(unitnum, '(a,2(I2,a),I2)')  "## Miller indexes: h = ", Miller_index(1), " ; k = ", Miller_index(2) &
            , " ; l = ", Miller_index(3)
        write(unitnum, '(a)') "## t | structure_factor(K,t)"
        do i =1, size(structure_factor)
            write(unitnum, format_style0) real(i,pr)*dt*conversion_factors(2), structure_factor(i)
        end do
    close(unitnum)

end subroutine write_structure_factor

subroutine write_observables(unitnum, time, energies, pressures, temperatures)
    real (pr), intent (in)              :: time, energies(2), pressures, temperatures
    integer(int_medium), intent (in)    :: unitnum

    write(unitnum, format_style0) time*conversion_factors(2), energies*conversion_factors(4), pressures*conversion_factors(6) &
        , temperatures

end subroutine write_observables

END MODULE writing2files