MODULE structureFactorMod
    use precisionMod
    use variablesMod
    use subroutinesMod
    implicit none

CONTAINS

subroutine get_reciprocal_vec(Miller_index, reciprocal_vec)
    integer, intent(in)     :: Miller_index(3)
    real(pr), intent(out)   :: reciprocal_vec(3)
    real(pr)                :: reciprocal_basis(3,3), reciprocal_latticeParam
    integer                 :: i

    select case (structure)
        case ("FCC")
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/-0.5_pr,  0.5_pr,  0.5_pr/)  ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/ 0.5_pr, -0.5_pr,  0.5_pr/)  ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/ 0.5_pr,  0.5_pr, -0.5_pr/)  ! reciprocal lattice basis vector 3
        case ("BCC")
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/0.0_pr, 0.5_pr, 0.5_pr/)      ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/0.5_pr, 0.0_pr, 0.5_pr/)      ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/0.5_pr, 0.5_pr, 0.0_pr/)      ! reciprocal lattice basis vector 3
        case ("random")
            print*, "Random structure selected -> Simple cubic basis set will be employed for structure factor calculation"
            reciprocal_latticeParam = 2._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/1._pr, 0._pr, 0._pr/)         ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/0._pr, 1._pr, 0._pr/)         ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/0._pr, 0._pr, 1._pr/)         ! reciprocal lattice basis vector 3
        case default
            reciprocal_latticeParam = 4._pr*pi/lattice_constant
            reciprocal_basis(:,1) = (/-0.5_pr,  0.5_pr,  0.5_pr/)   ! reciprocal lattice basis vector 1
            reciprocal_basis(:,2) = (/ 0.5_pr, -0.5_pr,  0.5_pr/)   ! reciprocal lattice basis vector 2
            reciprocal_basis(:,3) = (/ 0.5_pr,  0.5_pr, -0.5_pr/)   ! reciprocal lattice basis vector 3
    end select

    reciprocal_basis = reciprocal_basis * reciprocal_latticeParam

    reciprocal_vec = 0._pr
    do i =1, 3
        reciprocal_vec = reciprocal_vec + real(Miller_index(i),pr)*reciprocal_basis(:,i)
    end do

end subroutine get_reciprocal_vec

subroutine get_structure_factor(positions, structure_factor, reciprocal_vec)
    real(pr), intent(in)        :: positions(:,:), reciprocal_vec(3)
    real(pr), intent(inout)     :: structure_factor
    complex(pr)                 :: summation
    integer                     :: i

    summation = (0._pr,0._pr)

    do i = 1, num_atoms
        summation = summation + exp(CMPLX(0._pr,sum(reciprocal_vec*positions(:,i)), pr))
    end do

    structure_factor = abs(summation)*abs(summation)/(num_atoms*num_atoms)

end subroutine get_structure_factor

END MODULE structureFactorMod