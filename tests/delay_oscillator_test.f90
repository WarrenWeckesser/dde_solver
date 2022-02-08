!
! delay_oscillator_demo.f90
!
! Use DDE_SOLVER_M to solve the DDEs defined in delay_oscillator.f90.
!
! The program checks the value of the solution at t=5*pi/2, where the
! exact solution is 1.0.
!

PROGRAM delay_oscillator_test

USE DEFINE_delay_oscillator_DDEs
USE DDE_SOLVER_M

IMPLICIT NONE

INTEGER, DIMENSION(2) :: NVAR = (/NEQN, NLAGS/)

TYPE(DDE_SOL) :: SOL
TYPE(DDE_OPTS) :: OPTS

DOUBLE PRECISION, DIMENSION(1) :: LAGS
DOUBLE PRECISION, DIMENSION(1) :: p_
DOUBLE PRECISION, DIMENSION(2) :: TSPAN

DOUBLE PRECISION :: relerr, abserr, stoptime, finalerror, finaltol
DOUBLE PRECISION :: Pi

Pi = 3.1415926535897932385D0

! Intentionally too small to pass.
finaltol = 1.0D-21

! Set the parameters of the DDE
tau = Pi/2.0D0

! Set the solver parameters: relative error, abs. error, stop time
relerr = 1D-9
abserr = 1D-12
stoptime = 5*Pi/2
p_(1) = tau

! Initialize the array of lags
LAGS(1) = tau

TSPAN(1) = 0.0
TSPAN(2) = stoptime
OPTS = DDE_SET(RE=relerr, AE=abserr)

SOL = DDE_SOLVER(NVAR, delay_oscillator_ddes, LAGS, delay_oscillator_history, &
                 TSPAN, OPTIONS=OPTS)

finalerror = DABS(SOL%Y(SOL%NPTS, 1) - 1.0D0)

IF (finalerror > finaltol) THEN
    print *, 'finalerror = ', finalerror, ' exceeds desired max error of', finaltol
    ! NOTE: EXIT(status) is not standard Fortran!
    CALL EXIT(1)
END IF

END PROGRAM delay_oscillator_test
