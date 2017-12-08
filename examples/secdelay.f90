MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN):: T,Y,Z
    INTENT(OUT) :: DY
      DY(1) = Z(1,1)
    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y
    INTENT(IN)  :: T
    INTENT(OUT) :: Y
      Y(1) = T*T
    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE BETA(T,Y,BVAL)
  USE DDE_SOLVER_M
    DOUBLE PRECISION :: T,TVAL
    DOUBLE PRECISION, DIMENSION(:) :: Y
    DOUBLE PRECISION, DIMENSION(:) :: BVAL
    INTENT(IN):: T,Y
    INTENT(OUT) :: BVAL
    DOUBLE PRECISION, DIMENSION(NEQN) :: YOFT,DYOFT
      TVAL = T - T*T
      IF (TVAL <= 0D0) THEN
         CALL HISTORY(TVAL,YOFT)
!        YOFT(1) = TVAL * TVAL
      ELSE
!        USE THE SOLUTION QUEUE TO INTERPOLATE YOFT.
         CALL DDE_USER(TVAL,YOFT,DYOFT)
      END IF
      BVAL(1) = T - YOFT(1)
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE CHECK(T,Y,ERROR,ABSERR,RELERR)
    DOUBLE PRECISION :: T,X,FX,XGAMMA,FGAMMA,ABSERR,RELERR
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION :: EXACT,ERROR
    XGAMMA = 0.175487766624669D+01
    X = XGAMMA
    FGAMMA = X**9 / 9.0D0 - X**8 / 2.0D0 + (6.0D0/7.0D0) * X**7 -  &
             X**6 + X**5 - X**4 / 2.0D0 + X**3 / 3.0D0
    IF (T.LE.XGAMMA) THEN
       IF (T.LE.0.0D0) THEN
          EXACT = T * T
       ELSE
          EXACT = 0.0D0
       END IF
    ELSE
       X = T
       FX = X**9 / 9.0D0 - X**8 / 2.0D0 + (6.0D0 / 7.0D0) * X**7 - &
            X**6 + X**5 - X**4 / 2.0D0 + X**3 / 3.0D0
       EXACT = FX - FGAMMA
    ENDIF
!   COMPUTE THE ERROR OVERRUN.
    ERROR = ABS((Y(1)-EXACT)/(ABSERR+RELERR*ABS(EXACT)))
    RETURN
  END SUBROUTINE CHECK

END MODULE define_DDEs

!******************************************************************

PROGRAM secdelay 

! For t > 1, solve the DDE  y'(t) = y(t-t^2) for t >= 0,
!                           y(t) = t^2 for t <= 0.       
!
! The DDE is defined in the module define_DDEs. The problem is
! solved here with DKLAG6 and its output written to a file.
! Using the auxilary function getplotdata, it is then easy to
! import the data into Matlab and plot it. Interpolation is
! performed in subroutine BETA to approximate y(t-t^2).

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
  DOUBLE PRECISION :: ERROR,MAXERROR,ABSERR,RELERR,TVAL,TCHECK
  DOUBLE PRECISION, DIMENSION(NEQN) :: YCHECK
  DOUBLE PRECISION, DIMENSION(NEQN) :: YOFT,DYOFT

  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(12) :: FNAME='secdelay.dat'  

! Test illegal call:
  PRINT *,' Test with illegal call; should see an error message next.'
  TVAL = 0D0
  CALL DDE_USER(TVAL,YOFT,DYOFT)
  PRINT *,' '

  ABSERR = 1D-3
  RELERR = 1D-6
  OPTS = DDE_SET(RE=ABSERR,AE=RELERR)

  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,&
                  TSPAN=(/ 0D0, 5D0 /),OPTIONS=OPTS)

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

END PROGRAM secdelay
