module subrutinas
    use mtmod,    only: grnd    ! Mersenne -Twister RNG
    use precision
    use formats
    use mzranmod
    use funciones
    implicit none

    integer :: seed = -12345
    procedure(rng_interface), pointer :: get_rand => null()

    abstract interface
        function rng_interface() result(r)
            import :: pr
            real(kind=pr) :: r
        end function rng_interface
    end interface

contains

function ran2_wrapper() result(r)
    real(kind=pr) :: r
    r = ran2(seed)
end function ran2_wrapper

function rmzran_wrapper() result(r)
    real(kind=pr) :: r
    r = rmzran()
end function rmzran_wrapper

function grnd_wrapper() result(r)
    real(kind=pr) :: r
    r = grnd()
end function grnd_wrapper

subroutine set_rng(rng_name, user_seed)
    character(len=*), intent(in) :: rng_name
    integer, intent(in), optional :: user_seed

    if (present(user_seed)) seed = user_seed

    select case (trim(adjustl(rng_name)))
    case ("ran2")
        get_rand => ran2_wrapper
    case ("MZ")
        call mzran_init()
        get_rand => rmzran_wrapper
    case ("MT")
        get_rand => grnd_wrapper
    case default
        stop "Unknown RNG specified"
    end select
end subroutine set_rng


subroutine create_file_name(prefix, num, suffix, filename)
    character(len=8)                    :: fmt  ! Format descriptor
    character(len=12)                   :: x1   ! Temporary string for formatted real
    character(len=*)                    :: prefix, suffix
    character(len=:), allocatable       :: filename
    integer(kind=int_small), intent(in) :: num

    fmt = '(I4)'  ! Format integer
    write(x1, fmt) num  ! Convert real to string

    ! Trim spaces in formatted number
    x1 = adjustl(trim(x1))

    ! Concatenate strings
    filename = prefix // trim(x1) // suffix

end subroutine create_file_name

subroutine get_mean_square_displacement(N_iterations, dim, steps, mean_quadratic_displacement)
    integer(kind=int_large), intent(in)         :: N_iterations
    integer(kind=int_small), intent(in)         :: dim, steps
    real(kind=pr), intent(out)                  :: mean_quadratic_displacement
    real(kind=pr)                               :: rnd_num
    integer(kind=int_large), allocatable        :: lattice_position(:)
    integer(kind=int_large)                     :: i, j
    integer(kind=int_small)                     :: k

    allocate(lattice_position(dim))

    mean_quadratic_displacement = 0._pr
    do i = 1, N_iterations
        lattice_position = 0
        do j = 1, steps
            do k = 1, dim
                rnd_num = get_rand()
                if (rnd_num<0.5_pr) then
                    lattice_position(k) = lattice_position(k) + 1
                else if(rnd_num>0.5_pr) then
                    lattice_position(k) = lattice_position(k) - 1
                end if
            end do
        end do
        do k = 1, dim
            mean_quadratic_displacement = mean_quadratic_displacement + real(lattice_position(k)*lattice_position(k),pr)
        end do
    end do
    mean_quadratic_displacement =  mean_quadratic_displacement/real(N_iterations,pr)

    deallocate(lattice_position)

end subroutine get_mean_square_displacement

subroutine get_cuadrant_freq(N_iterations, dim, steps, cuadrants)
    integer(kind=int_large), intent(in)         :: N_iterations
    integer(kind=int_small), intent(in)         :: dim, steps
    real(kind=pr), allocatable, intent(out)     :: cuadrants(:)
    real(kind=pr)                               :: rnd_num
    integer(kind=int_large), allocatable        :: lattice_position(:)
    integer(kind=int_large)                     :: i, j
    integer(kind=int_small)                     :: k

    allocate(lattice_position(dim), cuadrants(2**dim))

    cuadrants = 0._pr
    do i = 1, N_iterations
        lattice_position = 0
        do j = 1, steps
            do k = 1, dim
                rnd_num = get_rand()
                if (rnd_num<0.5_pr) then
                    lattice_position(k) = lattice_position(k) + 1
                else if(rnd_num>0.5_pr) then
                    lattice_position(k) = lattice_position(k) - 1
                end if
            end do
        end do
        if (lattice_position(1)>0 .and. lattice_position(2)>=0 ) then
            cuadrants(1) = cuadrants(1) + 1._pr
        else if (lattice_position(1)>=0 .and. lattice_position(2)<0 ) then
            cuadrants(2) = cuadrants(2) + 1._pr
        else if (lattice_position(1)<=0 .and. lattice_position(2)>0 ) then
            cuadrants(3) = cuadrants(3) + 1._pr
        else if (lattice_position(1)<0 .and. lattice_position(2)<=0 ) then
            cuadrants(4) = cuadrants(4) + 1._pr
        end if

        ! select a random quadrant if the position is exactly at the origin
        if (lattice_position(1)==0 .and. lattice_position(2)==0) then
                rnd_num = get_rand()
                if (0._pr<=rnd_num .and. rnd_num<0.25_pr) then
                    cuadrants(1) = cuadrants(1) + 1._pr
                else if (0.25_pr<=rnd_num .and. rnd_num<0.5_pr) then
                    cuadrants(2) = cuadrants(2) + 1._pr
                else if (0.5_pr<=rnd_num .and. rnd_num<0.75_pr) then
                    cuadrants(3) = cuadrants(3) + 1._pr
                else if (0.75_pr<=rnd_num .and. rnd_num<1._pr) then
                    cuadrants(4) = cuadrants(4) + 1._pr
                end if
        end if
    end do
    cuadrants = cuadrants/N_iterations

    deallocate(lattice_position)

end subroutine get_cuadrant_freq

end module subrutinas
