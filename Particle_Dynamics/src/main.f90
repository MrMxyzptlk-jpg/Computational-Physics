
!***************************************************************
! Program: ex1.f90
! Purpose:
!
! Description:
!
! Input: All parameters are in ATOMIC UNITS. The program adimensionalized the problem in order to calculate, and re-dimnesionalizes the values when writing to files.
!
!
! Output:
!
!
!
! Room for improvement:  further modularization for comprehension and ease of use should be done in the 'subrutinas' modue. That is, to subdivide it into smaller modules in the 'modulos' folder. Also, the get_measurements() should be moved to observablesMod but requires too many modifications.
!
! Author: Jerónimo Noé Acito Pino
!***************************************************************
program ex1
    use precisionMod
    use phase_initializeMod
    use phase_runMod
    use phase_endMod
    implicit none

    real(pr)            :: CPU_t_start, CPU_t_end, CPU_elapsed_time

    call phase_initialize()

    CPU_t_start = omp_get_wtime()
    call phase_run()
    CPU_t_end = omp_get_wtime()
    CPU_elapsed_time = CPU_t_end - CPU_t_start

    call phase_end(CPU_elapsed_time)

end program ex1