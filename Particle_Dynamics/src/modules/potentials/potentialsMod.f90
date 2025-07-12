! Wrapper module with the different potentials
MODULE potentialsMod
    use Lennard_JonesMod
    use Coulomb_EwaldMod
    use Coulomb_ReactionFieldMod   ! Not fully implemented yet
    use potentialContributionMod
    implicit none

END MODULE potentialsMod