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

      DY(1) = Y(1) * Z(1,1) / T

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE BETA(T,Y,BVAL)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y
    DOUBLE PRECISION, DIMENSION(:) :: BVAL
    INTENT(IN):: T,Y
    INTENT(OUT) :: BVAL

      BVAL(1) = LOG(Y(1))

    RETURN
  END SUBROUTINE BETA

  SUBROUTINE CHECK(T,Y,ERROR)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION :: E1, E2, E3, EXACT, ERROR

      E1 = EXP(1D0)
      E2 = EXP(2D0)
      E3 = EXP(3D0 - EXP(1D0 - E1))
      IF (T <= E1) THEN
         EXACT = T
      ELSE
         IF (T <= E2) THEN
            EXACT = EXP(T/E1)
         ELSE
           IF (T <= E3) THEN
              EXACT = (E1/(3D0 - LOG(T)))**E1
           ELSE
              PRINT *, ' Solution not known for T > E3.'
              EXACT = Y(1)
           ENDIF
         ENDIF
      ENDIF
      ERROR = ABS((Y(1) - EXACT)/EXACT)

    RETURN
  END SUBROUTINE CHECK
END MODULE define_DDEs

!******************************************************************

PROGRAM exsd1 

! For t > 1, solve the DDE  y'(t) = y(t)*y(ln(y(t)))/t for t >= 1,
!                           y(t) = 1 for t <= 1.       
!
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
  DOUBLE PRECISION :: ERROR,MAXERROR,TCHECK
  DOUBLE PRECISION, DIMENSION(NEQN) :: YCHECK

  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(10) :: FNAME='exsd1.dat'  

  OPTS = DDE_SET(RE=1D-5,AE=1D-5)

  SOL = DDE_SOLVER(NVAR,DDES,BETA,(/ 1D0 /),&
                  (/ 1D0, 10D0 /),OPTIONS=OPTS)

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
        ! CALL CHECK(SOL%T(I),SOL%Y(I,1),ERROR)
        TCHECK = SOL%T(I)
        DO J = 1,NEQN
           YCHECK(J) = SOL%Y(I,J)
        END DO
        CALL CHECK(TCHECK,YCHECK,ERROR)      
        MAXERROR = MAX(ERROR,MAXERROR)
     END DO

     PRINT *,' '
     PRINT *,' Comparison to the analytical solution shows that'
     PRINT *,' at all steps the error is less than ',MAXERROR
     PRINT *,' '

  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  ENDIF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
   
  CALL RELEASE_ARRAYS(SOL,OPTS)
   
  STOP
END PROGRAM exsd1
