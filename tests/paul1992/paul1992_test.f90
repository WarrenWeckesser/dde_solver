!
! paul1992_demo.f90
!
! Fortran 90 program that will use DDE_SOLVER_M to solve the DDEs defined
! in the vector field: paul1992
!
! Reference:
! 
! C. A. H. Paul, Developing a delay differential equation solver,
! J. Appl. Num. Math., Vol. 9, pp. 403-414, 1992.
!

PROGRAM paul1992_demo

USE testing
USE DEFINE_paul1992_DDEs
USE DDE_SOLVER_M

IMPLICIT NONE

INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

TYPE(DDE_SOL) :: SOL
TYPE(DDE_OPTS) :: OPTS

DOUBLE PRECISION, DIMENSION(2) :: TSPAN
DOUBLE PRECISION :: e
DOUBLE PRECISION :: stoptime, expected

e = 2.718281828459045235360287D0

stoptime = 10.0D0

TSPAN(1) = 0.0D0
TSPAN(2) = stoptime
OPTS = DDE_SET(RE=1D-10, AE=1D-12)

SOL = DDE_SOLVER(NVAR, paul1992_ddes, paul1992_beta, paul1992_history, &
                 TSPAN, OPTIONS=OPTS)

expected = (e/(3.0D0 - DLOG(stoptime + 1.0D0))) ** e

IF (.not. is_close(SOL%Y(SOL%NPTS, 1), expected, 1D-9, 'Y(N)')) THEN
    ! NOTE: EXIT(status) is not standard Fortran!
    CALL EXIT(1)
END IF

END PROGRAM paul1992_demo
