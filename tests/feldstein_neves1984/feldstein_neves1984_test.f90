!
! feldstein_neves1984_test.f90
!
! Fortran 90 program that will use DDE_SOLVER_M to solve the DDEs defined
! in the vector field: feldstein_neves1984
!
! Reference:
! A. Feldstein and K. W. Neves, High order methods for state-dependent
! delay differential Equations with nonsmooth solutions, SIAM J. Numer. Anal.,
! Vol. 21, No. 5., pp. 844-863 (October, 1984)
!

program feldstein_neves1984_demo

use testing, only: is_close
use define_feldstein_neves1984_DDEs, only: &
        NEQN, NLAGS, &
        feldstein_neves1984_ddes, &
        feldstein_neves1984_history, &
        feldstein_neves1984_beta
use dde_solver_m, only: dde_solver, dde_set, dde_opts, dde_sol, release_arrays

implicit none

integer, dimension(2) :: nvar = (/NEQN, NLAGS/)

type(dde_sol) :: sol
type(dde_opts) :: opts

double precision, dimension(2) :: tspan
double precision :: stoptime, expected

stoptime = 2.0d0

tspan(1) = 0.0d0
tspan(2) = stoptime
opts = dde_set(re=1d-9, ae=1d-12)

sol = dde_solver(nvar, feldstein_neves1984_ddes, feldstein_neves1984_beta, &
                 feldstein_neves1984_history, tspan, options=opts)

if (sol%flag .ne. 0) then
    print *, 'dde_solver failed. sol%flag =', sol%flag
    call exit(1)
end if

! This expression for the expected value is only valid in the
! interval 1 <= t <= 2.
expected = (stoptime + 1.0d0)/4.0d0 + 0.5d0 + &
           (1.0d0 - dsqrt(2.0d0)/2.0d0)*dsqrt(stoptime + 1.0d0)

if (.not. is_close(sol%y(sol%npts, 1), expected, 1.0d-8, 'y(n)')) then
    ! NOTE: EXIT(status) is not standard Fortran!
    call exit(1)
end if

call release_arrays(sol, opts)

end program feldstein_neves1984_demo
