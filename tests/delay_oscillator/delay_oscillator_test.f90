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

use testing
use define_delay_oscillator_functions
use dde_solver_m

implicit none

! NEQN, NLAGS and NEF are from the define_delay_oscillator_functions module.
integer, dimension(3) :: nvar = (/NEQN, NLAGS, NEF/)

type(dde_sol) :: sol
type(dde_opts) :: opts

double precision, dimension(NLAGS) :: lags
double precision, dimension(2) :: tspan

double precision, parameter :: PI=3.1415926535897932385D0, PIby6=0.5235987755982988730771D0
double precision, dimension(3) :: expected_event_times = (/ PIby6, 5*PIby6, 2*PI + PIby6 /)
integer i
character(len=30) :: varname

! Initialize the lag to PI/2
lags(1) = PI/2.0D0

tspan(1) = 0.0
tspan(2) = 5.0D0*PI/2.0D0
opts = dde_set(re=1D-9, ae=1D-12)

sol = dde_solver(nvar, ddes, lags, history, &
                 tspan, event_fcn=event, options=opts)

if (.not. is_close(sol%y(sol%npts, 1), 1.0D0, 1.0D-9, 'y(N)')) then
    call exit(1)
end if

if (sol%ne .ne. 3) then
    print *, 'detected', sol%ne, ' events, expected 3'
    call exit(1)
end if

do i = 1, 3
    write(varname, "(A,I1,A)") "te(", i, ")"
    if (.not. is_close(sol%te(i), expected_event_times(i), 1.0D-9, varname)) then
        call exit(1)
    end if
    write(varname, "(A,I1,A)") "ye(", i, ")"
    if (.not. is_close(sol%ye(i, 1), 0.5D0, 1.0D-9, varname)) then
        call exit(1)
    end if
end do

end program delay_oscillator_test
