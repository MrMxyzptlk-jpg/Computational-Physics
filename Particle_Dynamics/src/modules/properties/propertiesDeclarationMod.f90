! Module to declare all properties of the particles and system
MODULE propertiesMod
    use precisionMod
    implicit none

    ! For particles' distinctions
    character(len=3), allocatable   :: symbols(:)

    ! For different integrators
    real(pr), dimension(:,:), allocatable   :: positions, Energies
    real(pr), dimension(:,:), allocatable   :: forces
    real(pr), dimension(:,:), allocatable   :: velocities
    real(pr), dimension(:,:), allocatable   :: previous_forces  ! Internal property for velocity-Verlet algorithm
    real(pr), dimension(:), allocatable     :: previous_E_potential     ! Internal property for Metropolis algorithm

    ! Observables
    real(pr), dimension(:), allocatable     :: Pressures, Temperatures
    real(pr), dimension(:), allocatable     :: pair_corr, meanSqrDisplacement, structure_factor

    ! For Coulomb interactions
    real(pr), allocatable   :: dipoles(:,:)
    real(pr), allocatable   :: charges(:)

END MODULE propertiesMod