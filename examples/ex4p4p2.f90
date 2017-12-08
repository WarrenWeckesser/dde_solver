MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1
  
CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY

    DY(1) = 2D0*Z(1,1)/(1D0 + Z(1,1)**9.65D0) - Y(1)
    
    RETURN
  END SUBROUTINE DDES

END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p4p2

! Example 4.4.2 in the SGT book.
!
! The DDE is defined in the module define_DDEs.  The problem
! is solved here with DDE_SOLVER and its output written to a
! file.  The auxilary function EX4P4P2.M imports the data into
! Matlab and plots it.


  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !
  !   NEQN = number of equations
  !   NLAGS = number of delays
  !
  ! are defined in the module define_DDEs as PARAMETERs so 
  ! they can be used for dimensioning arrays here. They are 
  ! passed to the solver in the array NVAR.

  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

  ! Constant delay
  DOUBLE PRECISION, PARAMETER :: LAG=2D0
  DOUBLE PRECISION, DIMENSION(NLAGS) :: DELAY=(/ LAG /)
  ! Constant history
  DOUBLE PRECISION, DIMENSION(NEQN) :: HISTORY= (/ 0.5D0 /)
  TYPE(DDE_SOL) :: SOL
  TYPE(DDE_OPTS) :: OPTS
  TYPE(DDE_INT) :: Y,YLAG

  DOUBLE PRECISION :: T0=0D0,TFINAL=100D0
  ! Prepare for output:
  INTEGER, PARAMETER :: NOUT=1000
  DOUBLE PRECISION, DIMENSION(NOUT) :: T

  ! Local variable:
  INTEGER :: I

  OPTS = DDE_SET(RE=1D-6,AE=1D-6,INTERPOLATION=.TRUE.)

  SOL = DDE_SOLVER(NVAR,DDES,DELAY,HISTORY, &
                   (/ T0,TFINAL /),OPTIONS=OPTS)

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN

    ! Form values for phase plane plot:
    T = (/ (LAG + (I-1)*((TFINAL - LAG)/(NOUT-1)), I=1,NOUT) /)
    Y = DDE_VAL(T,SOL)
    YLAG = DDE_VAL(T - LAG,SOL)
    ! Write the solution to a file for subsequent plotting
    ! in Matlab.
    OPEN(UNIT=6, FILE='ex4p4p2.dat')
    DO I = 1,NOUT
        WRITE(UNIT=6,FMT='(2D12.4)') Y%YT(I,1),YLAG%YT(I,1)
    ENDDO
      
    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'ex4p4p2.dat'."
    PRINT *,' '
    PRINT *,' These results can be accessed in Matlab and'
    PRINT *,' y(t) plotted against ylag = y(t - lag) by'
    PRINT *,' '
    PRINT *," >> [y,ylag] = ex4p4p2;"
    PRINT *,' '
  
  ELSE
  
    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG
  
  ENDIF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
  
  PRINT *, 'Call 1'
  CALL RELEASE_ARRAYS(SOL,OPTS)
  PRINT *, 'Call 2'
  CALL RELEASE_INT(Y)
  PRINT *, 'Call 3'
  CALL RELEASE_INT(YLAG)

  STOP
END PROGRAM ex4p4p2
