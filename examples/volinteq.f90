MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1
  DOUBLE PRECISION EPSTOL

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
    USE Q1DAMODULE
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,2*NLAGS) :: Z
    DOUBLE PRECISION A,B,EPS,ANS,ERR,RPAR
    INTEGER KF,IFLG,IPAR
    DIMENSION RPAR(1),IPAR(1)
    B = T*T*Y(1)
    A = T*Y(1)
    EPS = EPSTOL
    CALL Q1DA(F,A,B,EPS,ANS,ERR,KF,IFLG,RPAR,IPAR)
    DY(1) = (ANS - 1D0) / (T*T*T)
    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)
    DOUBLE PRECISION :: T
    !DOUBLE PRECISION, DIMENSION(NEQN,2) :: Y
    DOUBLE PRECISION, DIMENSION(2*NEQN) :: Y
    ! NEW: 06/28/2013:
    ! Return solution in 1,...,NEQN
    ! Derivative = NEQN+1,...,2*NEQN
    INTENT(IN)  :: T
    INTENT(OUT) :: Y
      !Y(1,1) = 1D0 / T
      !Y(1,2) = -1D0 / (T*T)
      Y(1) = 1D0 / T
      Y(2) = -1D0 / (T*T)
    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE BETA(T,Y,BVAL)
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION, DIMENSION(NLAGS) :: BVAL
      BVAL(1) = T
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE CHECK(T,Y,ERROR,ABSERR,RELERR)
    DOUBLE PRECISION :: T,ABSERR,RELERR
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION :: EXACT,ERROR
    EXACT = 1D0/T
!   COMPUTE THE ERROR OVERRUN.
    ERROR = ABS((Y(1)-EXACT)/(ABSERR+RELERR*ABS(EXACT)))
    WRITE(7,111) T,ERROR
    111 FORMAT(' T/Overrun = ',2D20.10)
    RETURN
  END SUBROUTINE CHECK

  SUBROUTINE F(S,RPAR,IPAR,FVAL)
!   INTEGRAND EVALUATION FOR Q1DA
    USE DDE_SOLVER_M
    IMPLICIT NONE
    DOUBLE PRECISION S,RPAR,FVAL
    INTEGER IPAR,IQUEUE
    DIMENSION RPAR(1),IPAR(1)
    DOUBLE PRECISION, DIMENSION(NEQN) :: YOFT,DYOFT
    !DOUBLE PRECISION, DIMENSION(NEQN,2) :: Y
    DOUBLE PRECISION, DIMENSION(2*NEQN) :: Y
!   EVALUATE THE INTEGRAND FOR Q1DA.
!   CHECK IF THE FIRST INTEGRATION STEP HAS BEEN COMPLETED.
!   IF NOT, USE THE HISTORY FUNCTION TO DEFINE Y AND Y'.
    CALL DDE_IQUEUE(IQUEUE)
    IF (IQUEUE == 0) THEN
       CALL HISTORY(S,Y)
       !FVAL = S*S*S*Y(1,1)*Y(1,2)
       FVAL = S*S*S*Y(1)*Y(2)
       RETURN
    ENDIF
    IF (S <= 1D0) THEN
!      FVAL = S*S*S*(1/S)*(-1/S**2)
       FVAL = -1D0
    ELSE
!      USE THE SOLUTION QUEUE TO INTERPOLATE YOFT.
       CALL DDE_USER(S,YOFT,DYOFT)
       FVAL = S*S*S*YOFT(1)*DYOFT(1)
    ENDIF
    RETURN
    END SUBROUTINE F
END MODULE define_DDEs

!******************************************************************

PROGRAM volinteq 

! Solves the integral equation
! y'(t) = [int{s^3 y(s) y'(s) ds, s=t*y(t),...,t*y^2(t)} - 1] / t^3
! y(t)  = 1
! The DDE is defined in the module define_DDEs. The problem is
! solved here with DKLAG6 and its output written to a file.
! Using the auxilary function getplotdata, it is then easy to
! import the data into Matlab and plot it.

  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !
  !   NEQN  = number of equations
  !   NLAGS = number of delays
  !
  ! are defined in the module define_DDEs as PARAMETERs so
  ! that they can be used for dimensioning arrays here.
  ! They are passed to the solver in the array NVAR.

  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

  TYPE(DDE_SOL) :: SOL ! Solution in Matlab style. 
  ! The fields of SOL are expressed in terms of the
  ! number of differential equations, NEQN, and the 
  ! number of output points, NPTS:

  !   SOL%NPTS         -- NPTS,number of output points.
  !
  !   SOL%T(NPTS)      -- values of independent variable, T.  
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T. 

  TYPE(DDE_OPTS) :: OPTS

  ! Local variables:
  INTEGER :: I,J
  DOUBLE PRECISION :: ERROR,MAXERROR,ABSERR,RELERR,TCHECK
  DOUBLE PRECISION, DIMENSION(NEQN) :: YCHECK
  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(13) :: FNAME='volinteq.dat'  
  CHARACTER(9) :: FNAMED='debug.ans'
  OPEN(UNIT=7, FILE=FNAMED)
 
  ABSERR = 1D-3
  RELERR = 1D-6
  EPSTOL = MIN(ABSERR,RELERR) / 100D0

  OPTS = DDE_SET(RE=ABSERR,AE=RELERR,NEUTRAL=.TRUE.)
  SOL  = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,&
                   TSPAN=(/ 1D0, 2D0 /),OPTIONS=OPTS)

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN
 
     ! Write the solution to a file for subsequent plotting
     ! in Matlab.
     OPEN(UNIT=6, FILE=FNAME)
     ! Set up format for exporting data.
     EXPORT = "(D12.4"//REPEAT(",D12.4",NEQN)//")"
     DO I = 1,SOL%NPTS
        WRITE(UNIT=6,FMT=EXPORT) SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
     ENDDO
      
     PRINT *,' Normal return from DDE_SOLVER with results'
     PRINT *," written to the file '",FNAME,"'."
     PRINT *,' '
     PRINT *,' These results can be accessed in Matlab'
     PRINT *,' from the folder holding this file by'
     PRINT *,' '
     PRINT *," >> [t,y] = getplotdata('",FNAME,"');"
     PRINT *,' '
     PRINT *,' and then plotted with, e.g., '
     PRINT *,' '
     PRINT *,' >> plot(t,y)'
     PRINT *,' '

     MAXERROR = 0D0
     DO I = 1,SOL%NPTS
        TCHECK = SOL%T(I)
        DO J =1,NEQN
           YCHECK(J) = SOL%Y(I,J)
        END DO 
        CALL CHECK(TCHECK,YCHECK,ERROR,ABSERR,RELERR)
        MAXERROR = MAX(ERROR,MAXERROR)
     ENDDO
     PRINT *,' '
     PRINT *,' Comparison to the analytical solution shows that'
     PRINT *,' at all steps the error overrun is less than ', &
     MAXERROR
     PRINT *,' '
     CALL PRINT_STATS(SOL)
  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  ENDIF

  STOP
END PROGRAM volinteq
