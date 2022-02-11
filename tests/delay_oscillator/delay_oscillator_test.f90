!
! delay_oscillator_demo.f90
!
! Use DDE_SOLVER_M to solve the DDEs defined in delay_oscillator.f90.
!
! The exact solution to the problem solved is y = sin(t).
! The program checks the value of the solution at t=5*pi/2, where the
! exact solution is 1.0.  It also checks that the three events where
! y = 0.5 in the interval 0 <= t <= 5*pi/2 are detected.
!

program delay_oscillator_test

use testing, only: is_close
use define_delay_oscillator_functions, only: ddes, history, event, NEQN, NLAGS, NEF
use dde_solver_m, only: dde_set, dde_solver, dde_sol, dde_opts, dde_int, &
                        dde_val, release_arrays, release_int

implicit none

! NEQN, NLAGS and NEF are from the define_delay_oscillator_functions module.
integer, dimension(3) :: nvar = (/NEQN, NLAGS, NEF/)

type(dde_sol) :: sol
type(dde_opts) :: opts
type(dde_int) :: yint

double precision, dimension(NLAGS) :: lags
double precision, dimension(2) :: tspan

double precision, parameter :: PI=3.1415926535897932385d0
double precision, dimension(3) :: tint, yint_expected
double precision, dimension(3) :: expected_event_times = &
    (/ PI/6.0d0, 5.0d0*PI/6.0d0, 2.0d0*PI + PI/6.0d0 /)
integer i
character(len=30) :: varname

! Initialize the lag to PI/2
lags(1) = PI/2.0d0

tspan(1) = 0.0
tspan(2) = 5.0d0*PI/2.0d0
opts = dde_set(re=1d-9, ae=1d-12, interpolation=.true.)

sol = dde_solver(nvar, ddes, lags, history, &
                 tspan, event_fcn=event, options=opts)

if (sol%flag .ne. 0) then
    print *, 'dde_solver failed. sol%flag =', sol%flag
    call exit(1)
end if

if (.not. is_close(sol%y(sol%npts, 1), 1.0d0, 1.0d-9, 'y(N)')) then
    call exit(1)
end if

if (sol%ne .ne. 3) then
    print *, 'detected', sol%ne, ' events, expected 3'
    call exit(1)
end if

do i = 1, 3
    write(varname, "(A,I1,A)") "te(", i, ")"
    if (.not. is_close(sol%te(i), expected_event_times(i), 1.0d-9, varname)) then
        call exit(1)
    end if
    write(varname, "(A,I1,A)") "ye(", i, ")"
    if (.not. is_close(sol%ye(i, 1), 0.5d0, 1.0d-9, varname)) then
        call exit(1)
    end if
end do

! Use the iterpolator to get the values at t=pi/3, t=pi/2 and t=3*pi/2.
! The expected values are sqrt(3)/2, 1 and -1, respectively.
tint = (/ PI/3.0d0, PI/2.0d0, 3.0d0*PI/2.0d0 /)
yint_expected = (/ dsqrt(3.0d0)/2.0d0, 1.0d0, -1.0d0 /)

yint = dde_val(tint, sol)

do i = 1, 3
    write(varname, "(A,I1,A)") "yint%yt(", i, ", 1)"
    if (.not. is_close(yint%yt(i, 1), yint_expected(i), 1.0d-9, varname)) then
        call exit(1)
    end if
end do

call release_arrays(sol, opts)
call release_int(yint)

end program delay_oscillator_test
