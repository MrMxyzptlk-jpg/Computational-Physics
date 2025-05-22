MODULE observables
    use precision
    use subrutinas
    implicit none

    CONTAINS

subroutine update_observables_absMagnetization(energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle &
    , magnetization_per_particle, binder)
    real(pr), intent(in)    :: Energy, magnetization
    real(pr), intent(inout) :: u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle, magnetization_per_particle, binder

    energy_per_particle = Energy/N_spins
    magnetization_per_particle = abs(magnetization)/N_spins

    u_avg = u_avg + energy_per_particle
    uSqr_avg = uSqr_avg + energy_per_particle*energy_per_particle

    m_avg = m_avg + magnetization_per_particle
    mSqr_avg = mSqr_avg + magnetization_per_particle*magnetization_per_particle

end subroutine update_observables_absMagnetization

subroutine update_observables_normalMagnetization(energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle &
, magnetization_per_particle, binder)
    real(pr), intent(in)    :: Energy, magnetization
    real(pr), intent(inout) :: u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle, magnetization_per_particle, binder

    magnetization_per_particle = magnetization/N_spins
    energy_per_particle = Energy/N_spins
    u_avg = u_avg + energy_per_particle
    uSqr_avg = uSqr_avg + energy_per_particle*energy_per_particle
    m_avg = m_avg + magnetization_per_particle
    mSqr_avg = mSqr_avg + magnetization_per_particle*magnetization_per_particle

end subroutine update_observables_normalMagnetization

subroutine update_observables_andBinder_normalMagnetization(energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg &
    , energy_per_particle, magnetization_per_particle, binder)
    real(pr), intent(in)    :: Energy, magnetization
    real(pr), intent(inout) :: u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle, magnetization_per_particle, binder
    real(pr)                :: m_sqr

    energy_per_particle = Energy/N_spins
    magnetization_per_particle = magnetization/N_spins

    u_avg = u_avg + energy_per_particle
    uSqr_avg = uSqr_avg + energy_per_particle*energy_per_particle

    m_sqr = magnetization_per_particle*magnetization_per_particle
    m_avg = m_avg + magnetization_per_particle
    mSqr_avg = mSqr_avg + m_sqr
    binder = binder + m_sqr*m_sqr

end subroutine update_observables_andBinder_normalMagnetization

subroutine update_observables_andBinder_absMagnetization(energy, magnetization, u_avg, uSqr_avg, m_avg, mSqr_avg &
    , energy_per_particle, magnetization_per_particle, binder)
    real(pr), intent(in)    :: Energy, magnetization
    real(pr), intent(inout) :: u_avg, uSqr_avg, m_avg, mSqr_avg, energy_per_particle, magnetization_per_particle, binder
    real(pr)                :: m_sqr

    energy_per_particle = Energy/N_spins
    magnetization_per_particle = abs(magnetization)/N_spins

    u_avg = u_avg + energy_per_particle
    uSqr_avg = uSqr_avg + energy_per_particle*energy_per_particle

    m_sqr = magnetization_per_particle*magnetization_per_particle
    m_avg = m_avg + magnetization_per_particle
    mSqr_avg = mSqr_avg + m_sqr
    binder = binder + m_sqr*m_sqr

end subroutine update_observables_andBinder_absMagnetization


END MODULE observables