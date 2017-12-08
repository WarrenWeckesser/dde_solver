MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=3,NLAGS=1,NEF=1

  ! Physical parameters
  DOUBLE PRECISION :: A1,A2,a,g,alpha,omegan,eta

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DYDT)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DYDT 

    ! Local variables
    DOUBLE PRECISION :: E,D,P,ylag1,ylag3
    E = Y(1)
    D = Y(2)
    P = Y(3)
    ylag1 = Z(1,1)       ! Retarded field.
    ylag3 = Z(3,1)       ! Retarded phase.

    IF (T < 200D0) THEN
       eta = 0D0
    ELSE
       eta = 0.2D0
    END IF

    DYDT(1) = (0.5D0)*(D + A2/(1D0 + a*E**2) - 1D0)*E &
                  + eta*ylag1*COS(ylag3 - P)

    DYDT(2) = g*(A1 - (1D0 + E**2)*D)

    DYDT(3) = omegan + (0.5D0*alpha)*(D + A2/(1D0 + a*E**2)) &
                       + eta*(ylag1/E)*SIN(ylag3 - P)

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE EF(T,Y,DYDT,Z,G)
  ! Inversion is max: dDdt=0. Not terminal. 
  ! dDdt is decreasing. 
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DYDt
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    DOUBLE PRECISION, DIMENSION(:) :: G
    INTENT(IN)  :: T,Y,DYDT,Z
    INTENT(OUT) :: G

    G(1) = A1 - (1D0 + Y(1)**2)*Y(2)

    RETURN
  END SUBROUTINE EF

END MODULE define_DDEs

!******************************************************************

PROGRAM laserex

  ! Self-Pulsing Semiconductor Laser with External Cavity
  
  ! Example provided by T.W. Carr.  Events are computed, 
  ! but not processed in this program.  The problem is 
  ! solved here with DDE_SOLVER and its output written to 
  ! a file.  The auxilary function LASEREX.M imports the 
  ! data into Matlab plots it in a manner equivalent to 
  ! Carr's solution with DDE23.M.
  !

  ! Variables and parameters in the model:
  !
  !  E = y(1)     Electric field magnitude.
  !  D = y(2)     Inversion.
  !  P = y(3)     Electric field phase.
  !  A1           Active pump.
  !  A2           Passive pump.
  !  a            Saturation ratio.
  !  g            Decay constant.
  !  alpha        Linewidth enhancement factor.
  !  omegan       External cavity phase.
  !  eta          Feedback strength.


  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !
  !   NEQN = number of equations
  !   NLAGS = number of delays
  !   NEF = number of event functions 
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

  TYPE(DDE_OPTS) :: OPTS

  DOUBLE PRECISION, DIMENSION(NEQN) :: HISTORY= &
  (/ 8.69863300579174D0, 0.11387116103411D0, 0D0 /)

  INTEGER :: I,J  ! Local variables

  DOUBLE PRECISION, PARAMETER :: DELAY=20D0,T0=0D0,TF=400D0
  ! For plotting y(t) against y(t - delay) in LASEREX.M, the
  ! mesh must have the following form. Note that the last point
  ! of the integration,TSPAN(NOUT), may not be exactly TF.
  DOUBLE PRECISION, PARAMETER :: REFINE=50,SPACING=DELAY/REFINE
  INTEGER, PARAMETER :: NOUT = 1 + (TF - T0)/SPACING
  DOUBLE PRECISION, DIMENSION(NOUT) :: &
  TSPAN = (/ (T0 + (I-1)*SPACING, I=1,NOUT) /)

  ! The physical parameters are global variables
  ! in define_DDEs.  Assign values here.
  A1 = 7D0
  A2 = -3.5D0  
  a = 2D0              
  g = 0.05D0          
  alpha = 5D0
  omegan = 0D0
  ! In the DDEs eta = 0 for T < 200 and eta = 0.2 for
  ! 200 <= T. Correspondingly, there is a jump at 200.

  ! There is an event function to locate the maximum of
  ! Inversion: dDdt=0. Not terminal. dDdt is decreasing. 
 
  OPTS = DDE_SET(RE=1D-5,AE=1D-9,DIRECTION=(/ -1 /),&
                 JUMPS=(/ 200D0 /))

  SOL = DDE_SOLVER(NVAR,DDES,(/ DELAY /),HISTORY,TSPAN,&
                   OPTIONS=OPTS,EVENT_FCN=EF) 

  ! Was the solver successful?
  IF(SOL%FLAG == 0) THEN

    PRINT *,' '    
    ! Print integration statistics:
    CALL PRINT_STATS(SOL)
 
    ! Write the solution to a file for subsequent 
    ! plotting in Matlab. 
    OPEN(UNIT=6, FILE='laserex.dat')
    DO I = 1,SOL%NPTS
      WRITE(UNIT=6,FMT='(4D15.6)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
    ENDDO
    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'laserex.dat'."
    PRINT *,' '
    PRINT *,' These results can be accessed in Matlab'
    PRINT *," and plotted as in T.W. Carr's program by"
    PRINT *,' '
    PRINT *," >> [t,y] = laserex;"
    PRINT *,' '
    PRINT *,' There were ',SOL%NE,' events.'
    PRINT *,' '
    PRINT *,' The times and corresponding solutions are:'
    PRINT *,' '
    DO I = 1,SOL%NE
      WRITE(*,FMT='(4D15.6)') SOL%TE(I),(SOL%YE(I,J),J=1,NEQN)
    ENDDO
    PRINT *,' '

  ELSE

    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG

  ENDIF
 
  CALL RELEASE_ARRAYS(SOL,OPTS)
   
  STOP
END PROGRAM laserex
