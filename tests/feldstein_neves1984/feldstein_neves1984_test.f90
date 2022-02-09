!
! feldstein_neves1984_demo.f90
!
! Fortran 90 program that will use DDE_SOLVER_M to solve the DDEs defined
! in the vector field: feldstein_neves1984
!
! Reference:
! A. Feldstein and K. W. Neves, High order methods for state-dependent
! delay differential Equations with nonsmooth solutions, SIAM J. Numer. Anal.,
! Vol. 21, No. 5., pp. 844-863 (October, 1984)
!

PROGRAM feldstein_neves1984_demo

USE DEFINE_feldstein_neves1984_DDEs
USE DDE_SOLVER_M

IMPLICIT NONE

INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

TYPE(DDE_SOL) :: SOL
TYPE(DDE_OPTS) :: OPTS

DOUBLE PRECISION, DIMENSION(2) :: TSPAN
DOUBLE PRECISION :: relerr, abserr, stoptime, expected, finalerror, finaltol

finaltol = 1.0D-8

! Set the solver parameters: relative error, abs. error, stop time
relerr = 1D-9
abserr = 1D-12
stoptime = 2.0D0

TSPAN(1) = 0.0D0
TSPAN(2) = stoptime
OPTS = DDE_SET(RE=relerr, AE=abserr)

SOL = DDE_SOLVER(NVAR, feldstein_neves1984_ddes, feldstein_neves1984_beta, &
                 feldstein_neves1984_history, TSPAN, OPTIONS=OPTS)

! This expression for the expected value is only valid in the
! interval 1 <= t <= 2.
expected = (stoptime + 1.0D0)/4.0D0 + 0.5D0 + &
           (1.0D0 - DSQRT(2.0D0)/2.0D0)*DSQRT(stoptime + 1.0D0)
finalerror = DABS(SOL%Y(SOL%NPTS, 1) - expected) / expected

IF (finalerror > finaltol) THEN
    print *, 'finalerror = ', finalerror, ' exceeds desired max error of', finaltol
    ! NOTE: EXIT(status) is not standard Fortran!
    CALL EXIT(1)
END IF

END PROGRAM feldstein_neves1984_demo
