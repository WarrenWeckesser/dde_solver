!
! paul1992_test.f90
!
! Fortran 90 program that will use DDE_SOLVER_M to solve the DDEs defined
! in the vector field: paul1992
!
! Reference:
! 
! C. A. H. Paul, Developing a delay differential equation solver,
! J. Appl. Num. Math., Vol. 9, pp. 403-414, 1992.
!

program paul1992_demo

use testing, only: is_close
use define_paul1992_ddes, only: &
        NEQN, NLAGS, &
        paul1992_ddes, &
        paul1992_history, &
        paul1992_beta
use dde_solver_m, only: dde_solver, dde_set, dde_opts, dde_sol, release_arrays

implicit none

double precision, parameter :: e=2.718281828459045235360287d0

integer, dimension(2) :: nvar = (/NEQN, NLAGS/)

type(dde_sol) :: sol
type(dde_opts) :: opts

double precision, dimension(2) :: tspan
double precision :: stoptime, expected

stoptime = 10.0d0

tspan(1) = 0.0d0
tspan(2) = stoptime
opts = dde_set(re=1d-10, ae=1d-12)

sol = dde_solver(nvar, paul1992_ddes, paul1992_beta, paul1992_history, &
                 tspan, options=opts)

if (sol%flag .ne. 0) then
    print *, 'dde_solver failed. sol%flag =', sol%flag
    call exit(1)
end if

expected = (e/(3.0d0 - dlog(stoptime + 1.0d0))) ** e

if (.not. is_close(sol%y(sol%npts, 1), expected, 1d-9, 'y(n)')) then
    ! note: exit(status) is not standard fortran!
    call exit(1)
end if

call release_arrays(sol, opts)

end program paul1992_demo
