! Wrapper module with the different potentials
MODULE potentialsMod
    use Lennard_JonesMod
    use Coulomb_EwaldMod
    use Coulomb_ReactionFieldMod   ! Not fully implemented yet
    use potentialContributionMod
    use potentialPointersMod
    implicit none

END MODULE potentialsMod