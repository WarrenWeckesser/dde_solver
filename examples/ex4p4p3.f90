MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=6,NLAGS=1,NEF=1
  ! Physical parameters:
  DOUBLE PRECISION, PARAMETER :: D=5D0,c=0.5D0,G=1D0,n=100 
  ! Flag that indicates whether to evaluate the ODE model 
  ! in subroutine DDES or the DDE model.
  LOGICAL :: useODEmodel

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
  ! x = y(1), y = y(2), lambda = y(3),
  ! I_1 = y(4), I_2 = y(5), I_3 = y(6).

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY

    IF (useODEmodel) THEN 
       ! Integrate the ODE model for T <= D:
       DY(1) = -Y(1)*Y(3) + G*Y(2)
       DY(2) = -DY(1)
       DY(4) = EXP(G*T)*Y(2)
       DY(5) = T*EXP(G*T)*Y(1)*Y(3)
       DY(6) = EXP(G*T)*Y(1)*Y(3)
       DY(3) = (c/n)*EXP(-G*T)*((DY(4)+DY(5))-G*(Y(4)+Y(5)))
    ELSE
       ! Integrate the DDE model for T >= D:
       DY(1) = -Y(1)*Y(3) + G*Y(2)
       DY(2) = -DY(1)
       DY(4) = EXP(G*T)*Y(2) - EXP(G*(T-D))*Z(2,1)
       DY(5) = D*EXP(G*T)*Y(1)*Y(3) - Y(6)
       DY(6) = EXP(G*T)*Y(1)*Y(3) - EXP(G*(T-D))*Z(1,1)*Z(3,1)
       DY(3) = (c/n)*EXP(-G*T)*((DY(4)+DY(5))-G*(Y(4)+Y(5)))
    ENDIF

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE EF(T,Y,DY,Z,G)
  ! Event function to recognize when to change from the 
  ! ODE model to the DDE model in subroutine DDES.

     DOUBLE PRECISION :: T
     DOUBLE PRECISION, DIMENSION(:) :: Y,DY
     DOUBLE PRECISION, DIMENSION(:,:) :: Z
     DOUBLE PRECISION, DIMENSION(:) :: G
     INTENT(IN) :: T,Y,DY,Z
     INTENT(OUT) :: G

     G(1) = T - D

     RETURN
  END SUBROUTINE EF

  SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
                  DIRECTION,ISTERMINAL,QUIT)
  ! Function to change a flag so that the DDE model will
  ! be evaluated in subroutine DDES instead of the ODE model.

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
       useODEmodel = .FALSE.
     END IF

    RETURN
  END SUBROUTINE CHNG

END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p4p3 

! Example 4.4.3 in the SGT book (multiple partnership model).

! The DDE is defined in the module define_DDEs.  The problem
! is solved here with DDE_SOLVER and its output written to a
! file.  The auxilary function EX4P4P3.M imports the data into
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

  ! Physical parameters D,c,G,n are defined in the module
  ! define_DDEs as PARAMETERS. 
  DOUBLE PRECISION, DIMENSION(NLAGS) :: LAGS=(/ D /)
  ! The simulation begins with an ODE model, so the 
  ! constant HISTORY is just y(0).
  DOUBLE PRECISION, DIMENSION(NEQN) :: Y0= &
    (/ 0.8D0*n,0.2D0*n,0D0,0D0,0D0,0D0 /)
  DOUBLE PRECISION, DIMENSION(2) :: TSPAN= (/0D0, 4D0*D /)

  ! Local variables:
  INTEGER :: I,J

  ! Initialize flag that indicates whether to evaluate 
  ! the ODE model in subroutine DDES or the DDE model.
  useODEmodel = .TRUE.

  SOL = DDE_SOLVER(NVAR,DDES,LAGS,Y0,TSPAN, & 
                   EVENT_FCN=EF,CHANGE_FCN=CHNG) 

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN

     ! Print integration statistics:
     CALL PRINT_STATS(SOL) 

     ! Write the first three components of the solution 
     ! to a file for subsequent plotting in Matlab.
     OPEN(UNIT=6, FILE='ex4p4p3.dat')
     DO I = 1,SOL%NPTS
        WRITE(UNIT=6,FMT='(4D12.4)') SOL%T(I),(SOL%Y(I,J),J=1,3)
     ENDDO
      
     PRINT *,' Normal return from DDE_SOLVER with results'
     PRINT *," written to the file 'ex4p4p3.dat'."
     PRINT *,' '
     PRINT *,' These results can be accessed in Matlab'
     PRINT *,' and plotted by'
     PRINT *,' '
     PRINT *," >> [t,y] = ex4p4p3;"
     PRINT *,' '
     PRINT *,' '

     PRINT *,' Kind of Event:'
     DO I = 1,SOL%NE
        IF(SOL%IE(I) == 1) THEN
           PRINT *,' Switched from ODEs to DDEs at t = ',SOL%TE(I)
        END IF
     END DO
     PRINT *,' '
 
  ELSE
 
     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG
 
  ENDIF

END PROGRAM ex4p4p3
