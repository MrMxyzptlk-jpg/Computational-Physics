
!***************************************************************
! Program: ex1.f90
! Purpose: simple simulation of particles interacting by Lennard-Jones or Coulomb potential.
!
! Description: Available 'integrators': velocity-Verlet, Monte-Carlo and Brownian. See further details in README.md
!
! Input: All parameters are in ATOMIC UNITS. The program adimensionalized the problem in order to calculate, and re-dimnesionalizes the values when writing to files.
!
!
! Output:
!
!
!
! Room for improvement: The Ewald implementation of the reciprocal space sum could be improved by having linked lists that index all good reciprocal vectors, thus reducing the sizes of all reciprocal arrays like eikx, eiky, eikz (after joining them into 1 array), kfac and reciprocal_charges.
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