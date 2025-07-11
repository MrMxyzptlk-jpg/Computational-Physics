! Wrapper module with the different updates subroutine for each integrator
MODULE updatePositionsMod
    use updatePositions_velVerletMod
    use updatePositions_MCMod
    use updatePositions_BDMod
    implicit none

END MODULE updatePositionsMod