MODULE define_DDEs

! Modified version of C2.

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=2,NLAGS=1,NEF=2
  INTEGER STATE1,STATE2

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
   
    DOUBLE PRECISION :: T,Y1,YLAG1
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
     
      Y1 = Y(1)
      YLAG1 = Z(1,1)
      DY(1) = -2D0*YLAG1
      IF (STATE2 .EQ. -1) THEN
         Y1 = -Y1
      END IF
      IF (STATE1 .EQ. -1) THEN
         YLAG1 = -YLAG1
      END IF
      DY(2) = (YLAG1 - Y1) / (1.0D0 + YLAG1)
    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y

     Y(1) = 1.0D0
     Y(2) = 0.5D0

    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE BETA(T,Y,BVAL)
    
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION, DIMENSION(NLAGS) :: BVAL

      BVAL(1) = T - Y(2)

    RETURN
  END SUBROUTINE BETA
  SUBROUTINE EF(T,Y,DY,Z,G)

    DOUBLE PRECISION :: T,YLAG1
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEF) :: G
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z

      YLAG1 = Z(1,1)
      G(1) = YLAG1
      G(2) = Y(1)
    
    RETURN
  END SUBROUTINE EF

 SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
     DIRECTION,ISTERMINAL,QUIT)

     INTEGER :: NEVENT
     LOGICAL, DIMENSION(NEF) :: ISTERMINAL
     INTEGER, DIMENSION(NEF) :: DIRECTION
     DOUBLE PRECISION :: TEVENT,HINIT
     DOUBLE PRECISION, DIMENSION(NEQN) :: YEVENT,DYEVENT
     LOGICAL :: QUIT

       IF (NEVENT == 1) THEN
          ! Reset the derivative switch for y_1(t).
          STATE1 = -STATE1
       ENDIF 

       IF (NEVENT == 2) THEN
          ! Reset the derivative switch for ylag_1(t).
          STATE2 = -STATE2
       ENDIF 

    RETURN
  END SUBROUTINE CHNG

END MODULE define_DDEs

!******************************************************************

PROGRAM c2g

! Program to solve a modified version of HENS C2. In this version
! root finding is used to locate the points at which y_1(t) = 0
! and ylag_1(t) = 0. Derivative switches are used for each to
! minimize the effect of the discontinuities caused by |y_1(t)|
! and |ylag_1(t)|.
!
! The DDE is defined in the module define_DDEs.  The 
! the problem is solved here with DKLAG6 and its output 
! written to a file.  Using the auxilary function 
! getplotdata, it is then easy to import the data into 
! Matlab and plot it.

  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !
  !   NEQN = number of equations
  !   NLAGS = number of delays
  !   NEF = number of event functions
  !
  ! are defined in the module define_DDEs as 
  ! PARAMETERs so that they can be used for 
  ! dimensioning arrays here. They are passed
  ! to the solver in the array NVAR.

  INTEGER, DIMENSION(3) :: NVAR = (/NEQN,NLAGS,NEF/)

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

  LOGICAL, DIMENSION(NEF) :: TERMINAL= (/ .FALSE.,.FALSE. /)
  INTEGER, DIMENSION(NEF) :: CROSS_AXIS=(/ 0,0 /)

  ! Local variables:
  INTEGER :: I,J

  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(10) :: FNAME='c2g.dat'  

  ! Format for writing events to an output file.
  CHARACTER(10+6*NEQN) :: EOUT
  ! EFNAME is the name of the output file.
  CHARACTER(13) :: EFNAME='c2gevents.dat'

! Initialize the derivative switches.
! Y1(T):
  STATE1 = 1
! YLAG1(T):
  STATE2 = 1

  OPTS = DDE_SET(ISTERMINAL=TERMINAL,DIRECTION=CROSS_AXIS,&
         RE=1D-5,AE=1D-5)

  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,&
        (/0D0,40D0/),EVENT_FCN=EF,OPTIONS=OPTS,CHANGE_FCN=CHNG)

  ! Was the solver successful?
  IF(SOL%FLAG == 0) THEN

    PRINT *,' '    
    ! Print integration statistics:
    CALL PRINT_STATS(SOL)
 
    ! Write the solution to a file for subsequent 
    ! plotting in Matlab.
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

     IF (SOL%NE > 0) THEN
        ! Display events to the screen and to a file for
        ! subsequent plotting in Matlab.
        OPEN(UNIT=8, FILE=EFNAME)
        EOUT = "(I5,D12.4"//REPEAT(",D12.4",NEQN)//")"
        PRINT *,' '
        PRINT *,' There were events that are displayed here and'
        PRINT *," written to the file '",EFNAME,"'."
        PRINT *,' '
        PRINT *,'  IE       TE                YE'
        DO I = 1,SOL%NE
            WRITE(*,EOUT) SOL%IE(I),SOL%TE(I),(SOL%YE(I,J),J=1,NEQN)
            WRITE(8,EOUT) SOL%IE(I),SOL%TE(I),(SOL%YE(I,J),J=1,NEQN)
        END DO
        PRINT *,' '
        PRINT *,' '
        PRINT *,' To access all the output in Matlab and plot the '
        PRINT *,' solution components as continuous curves with' 
        PRINT *,' events indicated by circles, use'
        PRINT *,' '
        PRINT *,' >> [t,y,te,ye,ie] = eventsc2g;'
        PRINT *,' '
     END IF
  ELSE

    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG

  ENDIF

  CALL RELEASE_ARRAYS(SOL,OPTS)
    
  STOP

END PROGRAM c2g

