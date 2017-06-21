MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1,NEF=3

  DOUBLE PRECISION :: R,MU,C

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T,YLAG
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
   
      YLAG = Z(1,1)

      IF (T <= 1D0-C) THEN
         DY(1) = -R*Y(1)*0.4D0*(1D0 - T)
      ELSEIF (T <= 1D0) THEN
         DY(1) = -R*Y(1)*( 0.4D0*(1D0 - T) + 10D0 - EXP(MU)*Y(1) )
      ELSEIF (T <= 2D0-C) THEN
         DY(1) = -R*Y(1)*( 10D0 - EXP(MU)*Y(1) )
      ELSE
         DY(1) = -R*EXP(MU)*Y(1)*( YLAG - Y(1) )
      END IF

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE EF(T,Y,DY,Z,G)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    DOUBLE PRECISION, DIMENSION(NEF) :: G

      G(1) = T - (1D0 - C)
      G(2) = T - 1D0
      G(3) = T - (2D0 - C)

    RETURN
  END SUBROUTINE EF

END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p7 

  ! Hoppensteadt and Waltman infection model, Exercise 4.7 in The 
  ! Book. Here event functions are used to force the solver to hit
  ! the three times at which the differential equation switches
  ! branches and the derivative is discontinuous at these times and
  ! to test the default values for isterminal and direction.

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
  ! that they can be used for dimensioning arrays here.
  ! They are passed to the solver in the array NVAR.

  INTEGER, DIMENSION(3) :: NVAR = (/NEQN,NLAGS,NEF/)

  TYPE(DDE_SOL) :: SOL
  ! Solution in Matlab style. 
  ! The fields of SOL are expressed in terms of the
  ! number of differential equations, NEQN, and the 
  ! number of output points, NPTS:

  !   SOL%NPTS         -- NPTS,number of output points.
  !
  !   SOL%T(NPTS)      -- values of independent variable, T.  
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T. 

  DOUBLE PRECISION, DIMENSION(NLAGS) :: LAGS=(/ 1D0 /)
  DOUBLE PRECISION, DIMENSION(NEQN) :: HISTORY=(/ 10D0 /)
  DOUBLE PRECISION, DIMENSION(2) :: TSPAN=(/ 0D0,10D0 /)

  ! Local variables:
  INTEGER :: I,J

  TYPE(DDE_OPTS) :: OPTS
  
  ! Format for writing solution to an output file.
  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(10) :: FNAME='export.dat'  

  ! Format for writing events to an output file.
  CHARACTER(10+6*NEQN) :: EOUT
  ! EFNAME is the name of the output file.
  CHARACTER(10) :: EFNAME='events.dat'  

  ! Define the physical parameters. They are passed
  ! as global variables via define_DDEs.
  R = 0.5D0
  MU = R/10D0
  C = 1D0/SQRT(2D0)

  OPTS = DDE_SET(RE=1D-5,AE=1D-8)

  SOL = DDE_SOLVER(NVAR,DDES,LAGS,HISTORY,TSPAN,&
                   EVENT_FCN=EF,OPTIONS=OPTS)

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

     IF(SOL%NE > 0) THEN
        ! Display events to the screen and to a file for
        ! subsequent plotting in Matlab.
        OPEN(UNIT=8, FILE=EFNAME)
        EOUT = "(I5,D12.4"//REPEAT(",D12.4",NEQN)//")"
        PRINT *,' '
        PRINT *,' There were events that are displayed here and'
        PRINT *," written to the file '",EFNAME,"'."
        PRINT *,' '
        PRINT *,'  IE       TE         YE'
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
        PRINT *,' >> [t,y,te,ye,ie] = ex4p7;'
        PRINT *,' '

     END IF

  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  ENDIF
  
  STOP
END PROGRAM ex4p7
