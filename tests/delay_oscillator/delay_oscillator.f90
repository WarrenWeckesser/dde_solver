!
! delay_oscillator.f90
!

module define_delay_oscillator_functions

    implicit none
    integer, parameter :: NEQN=1, NLAGS=1, NEF=1

contains

    subroutine ddes(t, x, Zlags, vf)
    ! Arguments
    double precision, intent(in) :: t
    double precision, dimension(:), intent(in) :: x
    double precision, dimension(:,:), intent(in) :: Zlags
    double precision, dimension(:), intent(out) :: vf

    vf(1) = -Zlags(1, 1)

    return
    end subroutine ddes

    subroutine history(t, x)
    double precision, intent(in) :: t
    double precision, dimension(:), intent(out) :: x

    x(1) = sin(t)

    return
    end subroutine history

    subroutine event(t, y, dy, z, g)
    ! Arguments
    double precision, intent(in) :: t
    double precision, dimension(:), intent(in) :: y, dy
    double precision, dimension(:,:), intent(in) :: z
    double precision, dimension(:), intent(out) :: g

    g(1) = y(1) - 0.5D0

    return
    end subroutine event

end module define_delay_oscillator_functions
