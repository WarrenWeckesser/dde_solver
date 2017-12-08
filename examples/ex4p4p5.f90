MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=2,NLAGS=1,NEF=2

  INTEGER :: STATE

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION :: YLAG
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN):: T,Y,Z
    INTENT(OUT) :: DY
   
 !  Physical parameters
    DOUBLE PRECISION, PARAMETER :: gamma=0.248D0,beta=1D0, &
                                   A=0.75D0,omega=1.37D0

    YLAG = Z(1,1)

    DY(1) = Y(2)       
    DY(2) = SIN(Y(1)) - STATE*gamma*COS(Y(1)) - beta*YLAG &
              + A*SIN(omega*T + ASIN(gamma/A))

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE EF(T,Y,DY,Z,G)

     DOUBLE PRECISION :: T
     DOUBLE PRECISION, DIMENSION(:) :: Y,DY
     DOUBLE PRECISION, DIMENSION(:,:) :: Z
     DOUBLE PRECISION, DIMENSION(:) :: G
     INTENT(IN) :: T,Y,DY,Z
     INTENT(OUT) :: G

     G = (/ Y(1), ABS(Y(1)) - ASIN(1D0) /)

     RETURN
  END SUBROUTINE EF

  SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
     DIRECTION,ISTERMINAL,QUIT)

     INTEGER :: NEVENT
     INTEGER, DIMENSION(:) :: DIRECTION
     DOUBLE PRECISION :: TEVENT,HINIT
     DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
     LOGICAL :: QUIT
     LOGICAL, DIMENSION(:) :: ISTERMINAL
     INTENT(IN) :: NEVENT,TEVENT
     INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,&
                      ISTERMINAL,QUIT

       IF (NEVENT == 1) THEN
          ! Restart the integration with initial values
          ! that correspond to a bounce of the suitcase.
          STATE = -STATE
          YEVENT(1) = 0.0D0
          YEVENT(2) = 0.913*YEVENT(2)
          DIRECTION(1) = - DIRECTION(1)
     ! ELSE
     !    The suitcase fell over, NEVENT = 2. The integration
     !    could be terminated by QUIT = .TRUE., but this 
     !    event is already a terminal event.
       ENDIF 

    RETURN
  END SUBROUTINE CHNG

END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p4p5 

! Example 4.4.4 in the SGT book.

! The DDE is defined in the module define_DDEs.  The problem
! is solved here with DDE_SOLVER and its output written to a
! file.  The auxilary function EX4P4P5.M imports the data into
! Matlab and plots it.

  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !
  !   NEQN  = number of equations
  !   NLAGS = number of delays
  !   NEF   = number of event functions 
  !
  ! are defined in the module define_DDEs as PARAMETERs so 
  ! they can be used for dimensioning arrays here. They are 
  ! passed to the solver in the array NVAR.

  INTEGER, DIMENSION(3) :: NVAR = (/NEQN,NLAGS,NEF/)

  TYPE(DDE_SOL) :: SOL 
  ! The fields of SOL are expressed in terms of the
  ! number of differential equations, NEQN, and the 
  ! number of output points, NPTS:

  !   SOL%NPTS         -- NPTS,number of output points.
  !
  !   SOL%T(NPTS)      -- values of independent variable, T.  
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T. 
  
  ! When there is an event function, there are fields

  !   SOL%NE           -- NE, number of events.
  !   SOL%TE(NE)       -- locations of events
  !   SOL%YE(NE,NEQN)  -- values of solution at events
  !   SOL%IE(NE)       -- identifies which event occurred

  TYPE(DDE_OPTS) :: OPTS

  ! Local variables:
  INTEGER :: I,J

  ! Prepare output points
  INTEGER, PARAMETER :: NOUT=1000
  DOUBLE PRECISION, PARAMETER :: T0=0D0,TFINAL=12D0
  DOUBLE PRECISION, DIMENSION(NOUT) :: TSPAN= &
  (/ (T0+(I-1)*((TFINAL - T0)/(NOUT-1)), I=1,NOUT) /)

  ! Initialize global variable that governs the
  ! form of the DDEs.
  STATE = 1

  OPTS = DDE_SET(RE=1D-5,DIRECTION=(/ -1,0 /), &
                 ISTERMINAL=(/ .FALSE.,.TRUE. /) )

  SOL = DDE_SOLVER(NVAR,DDES,(/ 0.1D0 /),(/ 0D0,0D0 /), &
         TSPAN,OPTIONS=OPTS,EVENT_FCN=EF,CHANGE_FCN=CHNG) 

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN
 
     ! Write the solution to a file for subsequent plotting
     ! in Matlab.
     OPEN(UNIT=6, FILE='ex4p4p5.dat')
     DO I = 1,SOL%NPTS
        WRITE(UNIT=6,FMT='(3D12.4)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
     ENDDO
      
     PRINT *,' Normal return from DDE_SOLVER with results'
     PRINT *," written to the file 'ex4p4p5.dat'."
     PRINT *,' '
     PRINT *,' These results can be accessed in Matlab'
     PRINT *,' and plotted in a phase plane by'
     PRINT *,' '
     PRINT *," >> [t,y] = ex4p4p5;"
     PRINT *,' '
     PRINT *,' '

     PRINT *,' Kind of Event:'
     DO I = 1,SOL%NE
        IF(SOL%IE(I) == 1) THEN
           PRINT *,' A wheel hit the ground at',SOL%TE(I)
        ELSE
           PRINT *,' The suitcase fell over at',SOL%TE(I)
        END IF
     END DO
     PRINT *,' '
 
  ELSE
 
     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG
 
  ENDIF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
  
  CALL RELEASE_ARRAYS(SOL,OPTS)
    
  STOP
END PROGRAM ex4p4p5
