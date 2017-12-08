MODULE define_DDEs

  IMPLICIT NONE

  INTEGER, PARAMETER :: NEQN=4,NLAGS=2

  ! Physical parameters assigned in the main program.
  DOUBLE PRECISION :: tau,omega

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY
   
    ! Physical parameters
    DOUBLE PRECISION :: A=0.330D0,d=0.006D0,lambda=0.308D0,&
                        gamma=0.040D0,epsilon=0.060D0

    ! Local variables
    DOUBLE PRECISION :: S,E,I,R,Noft,Itau,Somega,Eomega,&
    Iomega,Romega,Nomega,dSdt,dEdt,dIdt,dRdt

    ! Local solution variables and delayed solution values:
    S = Y(1)
    E = Y(2)
    I = Y(3)
    R = Y(4)
    Itau   = Z(3,1)
    Somega = Z(1,2)
    Eomega = Z(2,2)
    Iomega = Z(3,2)
    Romega = Z(4,2)

    Noft = S + E + I + R
    Nomega = Somega + Eomega + Iomega + Romega

    ! Local derivatives:
    dSdt = A - d*S - lambda*((S*I)/Noft) + gamma*Itau*exp(-d*tau)
    dEdt = lambda*((S*I)/Noft) -                                 &
           lambda*((Somega*Iomega)/Nomega)* exp(-d*omega) - d*E
    dIdt = lambda*((Somega*Iomega)/Nomega)*exp(-d*omega) -       &
           (gamma+epsilon+d)*I
    dRdt = gamma*I - gamma*Itau*exp(-d*tau) - d*R

    ! Derivatives for the integrator:
    DY = (/ dSdt, dEdt, dIdt, dRdt /)

    RETURN
  END SUBROUTINE DDES

END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p4p1 

! Example 4.4.1 in the SGT book.
!
! The DDE is defined in the module define_DDEs.  The problem
! is solved here with DDE_SOLVER and its output written to a
! file.  The auxilary function EX4P4P1.M imports the data into
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

  ! Constant delays and history:
  DOUBLE PRECISION, DIMENSION(NLAGS) :: DELAYS=&
  (/ 42D0, 0.15D0 /)
  DOUBLE PRECISION, DIMENSION(NEQN) :: HISTORY=&
  (/ 15D0, 0D0, 2D0, 3D0 /)

  TYPE(DDE_SOL) :: SOL
  !
  !   SOL%NPTS         -- NPTS,number of output points.
  !
  !   SOL%T(NPTS)      -- values of independent variable, T.  
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T. 

  INTEGER :: I,J ! Local variables

  ! These physical parameters are global variables
  ! in define_DDEs.  Assign their values here.
  tau = 42D0
  omega = 0.15D0

  SOL = DDE_SOLVER(NVAR,DDES,DELAYS,HISTORY,TSPAN=(/ 0D0,350D0 /))

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN

     ! Write the solution to a file for subsequent plotting
     ! in Matlab.
     OPEN(UNIT=6, FILE='ex4p4p1.dat')
     DO I = 1,SOL%NPTS
        WRITE(UNIT=6,FMT='(5D12.4)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
     ENDDO
      
     PRINT *,' Normal return from DDE_SOLVER with results'
     PRINT *," written to the file 'ex4p4p1.dat'."
     PRINT *,' '
     PRINT *,' These results can be accessed in Matlab'
     PRINT *,' and plotted by'
     PRINT *,' '
     PRINT *," >> [t,y] = ex4p4p1;"
     PRINT *,' '

  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  ENDIF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)

END PROGRAM ex4p4p1
