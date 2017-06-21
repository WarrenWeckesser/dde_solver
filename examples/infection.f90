MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=3,NLAGS=3,NEF=1
! Global variables: 
  DOUBLE PRECISION :: M,SIGMA,S0,R0,T0
  INTEGER :: STATE,EXAMPLE

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
!   Define the derivatives for the DDEs.
    USE DDE_SOLVER_M
    DOUBLE PRECISION :: T,TAU,TVAL,STAU1,STAU2,STAU3,STAU4
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
    DOUBLE PRECISION, DIMENSION(NEQN) :: YOFT,DYOFT
!   Local variables:
!   Y(1) = tau(t)
!   Y(2) = S(t)
!   Y(3) = integral(rho(s)*I0(s),s=0...t)
    IF (STATE == 1) THEN
!      tau'(t) = 0
       DY(1) = 0D0
!      S'(t) = -r(t) * I(t) * S(t):
       DY(2) = -RFUN(T) * I0FUN(T) * Y(2)
!      (Accumulated dosage)' = rho(t)*I0(t):
       DY(3) = RHOFUN(T) * I0FUN(T)
    ELSE
!      Z contains the delayed solution values. Its contents
!      are:
!        tau(t-sigma)      tau(tau(t))    tau(tau(t)-sigma)
!        S(tau-sigma)      S(tau(t))      S(tau(t)-sigma)
!        y_3(tau-sigma)    y_3(tau(t))    y_3(tau(t)-sigma)
!      The definition of tau'(t) requires both I(t) and
!      I(tau(t)). We need to approximate the following
!      quantities using subroutine DDE_USER.
!        S1 = S(tau(t-sigma))
!        S2 = S(tau(t))
!        S3 = S(tau(tau(t)-sigma))
!        S4 = S(tau(tau(t)))
       TVAL = Z(1,1)                      ! tau(t-sigma)
       IF (TVAL < 0D0) THEN
          STAU1 = S0
       ELSE
          CALL DDE_USER(TVAL,YOFT,DYOFT)
          STAU1 = YOFT(2)                 ! S(tau(t-sigma))
       ENDIF
       STAU2 = Z(2,2)                     ! S(tau(t))
       TVAL = Z(1,3)                      ! tau(tau(t)-sigma)
       IF (TVAL < 0D0) THEN
          STAU3 = S0
       ELSE
          CALL DDE_USER(TVAL,YOFT,DYOFT)
          STAU3 = YOFT(2)                 ! S(tau(tau(t)-sigma))
       ENDIF
       TAU = Y(1)                         ! tau(t)
       TVAL = Z(1,2)                      ! tau(tau(t))
       IF (TVAL < 0D0) THEN
          STAU4 = S0
       ELSE
          CALL DDE_USER(TVAL,YOFT,DYOFT)
          STAU4 = YOFT(2)                 ! S(tau(tau(t)))
       ENDIF
!      tau'(t) = (rho(t) * I(t)) / (rho(tau(t)) * I(tau(t))):
       DY(1) = (RHOFUN(T) * IFUN(T,STAU2,STAU1)) /   &
               (RHOFUN(TAU) * IFUN(TAU,STAU4,STAU3))
!      S'(t) = -r(t) * I(t) * S(t):
       DY(2) = -RFUN(T) * IFUN(T,STAU2,STAU1) * Y(2)
!      (Accumulated dosage)' = rho(t)*I0(t):
       DY(3) = RHOFUN(T) * I0FUN(T)
    ENDIF
    RETURN
  END SUBROUTINE DDES

  SUBROUTINE BETA(T,Y,BVAL)
!   Define the delay times.
    USE DDE_SOLVER_M
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION, DIMENSION(NLAGS) :: BVAL
!     The delay times are t-sigma, tau(t), and tau(t)-sigma:
!     t-sigma:
      BVAL(1) = T - SIGMA
!     tau(t):
      BVAL(2) = Y(1)
!     tau(t)-sigma:
      BVAL(3) = Y(1) - SIGMA
    RETURN
  END SUBROUTINE BETA

  SUBROUTINE EF(T,Y,DY,Z,G)
!     If STATE = 1, the minimum threshold time t_0
!     corresponds to the event G(1) = 0. If STATE = 2,
!     the time t_0 + sigma corresponds to the event
!     G(1) = 0.
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y,DY
    DOUBLE PRECISION, DIMENSION(NEF) :: G
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: Z
      IF (STATE == 1) THEN
         G(1) = Y(3) - M 
      ELSE
         G(1) = T - (T0 + SIGMA)
      ENDIF
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
          IF (STATE == 1) THEN
             STATE = 2
             T0 = TEVENT
             PRINT *,' Computed t0   = ', TEVENT
             IF (EXAMPLE == 1) &
             PRINT *,' Analytical t0 = ', (2D0-SQRT(2D0))/2D0 
             IF (EXAMPLE == 2) &
             PRINT *,' Approximate check t0 = ', 0.357403D0
          ELSE IF (STATE == 2) THEN
             PRINT *,' Computed t_0 + sigma = ', TEVENT
             STATE = 3
          ENDIF 
       ENDIF 
    RETURN
  END SUBROUTINE CHNG

 FUNCTION I0FUN(T) RESULT(I0)
!   Define I_0(t).
    DOUBLE PRECISION :: T,I0
    IF (T <= -1D0) THEN
       I0 = 0D0
    ELSEIF (T <= 0D0) THEN
       I0 = 0.4D0*(1D0 + T)
    ELSEIF (T <= 1D0) THEN
       I0 = 0.4D0*(1D0 - T)
    ELSE
       I0 = 0D0
    ENDIF
 END FUNCTION I0FUN

 FUNCTION RHOFUN(T) RESULT(RHO)
!   Define rho(t).
    DOUBLE PRECISION T,RHO
      RHO = 1D0
      IF (EXAMPLE == 2) RHO = EXP(-T)
 END FUNCTION RHOFUN

 FUNCTION RFUN(T) RESULT(R)
!   Define r(t).
    DOUBLE PRECISION :: T,R
      R = R0
      IF (EXAMPLE == 2) R = R0 * (1D0 + SIN(5D0*T))
 END FUNCTION RFUN

 FUNCTION IFUN(T,S1,S2) RESULT(I)
!   Define I(t).
    DOUBLE PRECISION :: T,I,S1,S2
    IF (STATE == 1) THEN
!      I(t) = I_0(t):
       I = I0FUN(T)
    ELSEIF (STATE == 2) THEN
!      I(t) = I_0(t) + S_0 - S(tau(t)):
       I = I0FUN(T) + S0 - S1
    ELSE
!      I(t) = S(tau(t-sigma)) - S(tau(t)):
       I = S2 - S1
    ENDIF
 END FUNCTION IFUN

END MODULE define_DDEs

!******************************************************************

PROGRAM infection 

! This program solves Example 1-2 from "Numerical Solution of a
! problem in the theory of epidemics" by F.C. Hoppensteadt and
! Z. Jackiewicz. An event function is used to locate the
! mimimum threshold time T0 and to locate the time T0+SIGMA.
! Details of how this is done are found in the document
! infection.pdf.

! The DDEs are defined in the module define_DDEs. The problem is
! solved here with DDE_SOLVER for four values of a parameter R0
! and the output written to files. The auxilary Matlab function 
! infection.m plots the solution S(t) and makes the data available
! as matrices for study in Matlab. The auxiliary Matlab function
! infectionIoft.m plots the function I(t).

  USE define_DDEs
  USE DDE_SOLVER_M

  IMPLICIT NONE

  ! The quantities
  !   NEQN  = number of equations
  !   NLAGS = number of delays
  ! are defined in the module define_DDEs as PARAMETERs so
  ! that they can be used for dimensioning arrays here.
  ! They are passed to the solver in the array NVAR.
  INTEGER, DIMENSION(3) :: NVAR = (/NEQN,NLAGS,NEF/)

  TYPE(DDE_SOL) :: SOL ! Solution structure 
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
  TYPE(DDE_INT) :: IPLOT

  ! Local variables:
  INTEGER :: I,J,RVAL
  ! Variables for plotting I(t):
  INTEGER, PARAMETER :: NPLOT = 101
  DOUBLE PRECISION :: IT
  DOUBLE PRECISION, DIMENSION(NPLOT) :: T,TAU,S2,T1,S1

  OPEN(UNIT=6, FILE='infection.dat')
  WRITE(UNIT=6,FMT="(I12)") NEQN
  OPEN(UNIT=7, FILE='Ioft.dat')

  ! Solve for four values of R0:
  DO RVAL = 1,4
    ! Define the problem parameters using global variables:
    R0 = 0.1D0*(RVAL+1)
    PRINT *,' '
    PRINT *,' Case R0 = ', R0
    M = 0.1D0
    SIGMA = 1D0
    S0 = 10D0
!   To solve Example 2 in the paper change the next line to
!   EXAMPLE = 2.
    EXAMPLE = 1
    IF (EXAMPLE == 1) PRINT *, ' Results for Example 1 follow.'
    IF (EXAMPLE == 2) PRINT *, ' Results for Example 2 follow.'
    
    ! The switch, STATE, is used to distinguish the
    ! the different phases of the integration. It
    ! is initialized to 1. It is changed to 2 when
    ! the minimum threshold time T0 is located. It
    ! is changed to 3 at time T0+SIGMA.

    STATE = 1
    
    ! An unusual aspect of this problem is that we cannot
    ! automatically track discontinuities because one of the
    ! delays, namely tau(t), is identically 0 on [0,T0].
    ! The discontinuities of interest occur at T0 and at
    ! T0+SIGMA; both are handled using events functions.

    ! Note that interpolation is to be done and stringent 
    ! tolerances are specified. However, default error
    ! tolerances yields an accurate solution. Interpolation
    ! is done to illustrate how to plot auxiliary functions,
    ! in this case I(t).

    OPTS = DDE_SET(TRACK_DISCONTINUITIES=.FALSE.,&
                   ISTERMINAL=(/.FALSE./),       &
                   DIRECTION=(/+1/),             &
                   INTERPOLATION=.TRUE.,         &
                   RE=1D-8,AE=1D-8)
    
    ! Note the constant history function is supplied as a
    ! vector.
    SOL  = DDE_SOLVER(NVAR,DDES,BETA,(/ 0D0,S0,0D0 /),  &
                      (/ 0D0,8D0 /),OPTIONS=OPTS,       &
                      EVENT_FCN=EF,CHANGE_FCN=CHNG)

    ! If the solver was successful, write the data necessary
    ! to plot S(t) and I(t) to files for use by the auxiliary
    ! Matlab functions infection.m and infectionIoft.m.

    IF (SOL%FLAG == 0) THEN
       WRITE(UNIT=6,FMT="(I12)") SOL%NPTS
       DO I = 1,SOL%NPTS
         WRITE(UNIT=6,FMT="(D12.4)") SOL%T(I),&
                           (SOL%Y(I,J),J=1,NEQN)
       END DO
       ! Generate the data to plot I(t).
       T = (/ ((I-1)*(8D0/(NPLOT-1)), I=1,NPLOT) /)
       IPLOT = DDE_VAL(T, SOL)
       TAU = IPLOT%YT(:,1)
       TAU = MAX(TAU,0D0)                  ! tau(t)
       IPLOT = DDE_VAL(TAU, SOL)
       S2 = IPLOT%YT(:,2)              ! S(tau(t))
       T1 = T - SIGMA               ! t-sigma
       T1 = MAX(T1,0D0)
       IPLOT = DDE_VAL(T1, SOL)
       S1 = IPLOT%YT(:,1)        ! tau(t-sigma)
       S1 = MAX(S1,0D0)
       IPLOT = DDE_VAL(S1, SOL)
       S1 = IPLOT%YT(:,2)        ! S(tau(t)-sigma)
       ! As infectionIoft.m is coded, a second dummy output
       ! variable is needed here.
       WRITE(UNIT=7,FMT="(2I12)") NPLOT,0
       DO I = 1,NPLOT
          IF (T(I) >= T0+SIGMA) THEN
             ! I(t) = S(tau(t-sigma)) - S(tau(t)):
             IT = S1(I) - S2(I)
          ELSEIF (T(I) >= T0) THEN
             ! I(t) = I0(t) + S0 - S(tau(t)):
             IT = I0FUN(T(I)) + S0 - S2(I)
          ELSE
             ! I(t) = I0(t):
             IT = I0FUN(T(I))
          ENDIF
          IT = MAX(IT,0D0)
          WRITE(UNIT=7,FMT="(2D12.4)") T(I),IT
       END DO
    ELSE
      ! If the solver was not successful, terminate with
      ! an error message.
      PRINT *,' Abnormal return with FLAG = ',SOL%FLAG
      STOP
    END IF

    ! Integration statistics:
    CALL PRINT_STATS(SOL)
   
!   Release the solution and interpolation structure arrays.
    CALL RELEASE_ARRAYS(SOL,OPTS)
    CALL RELEASE_INT(IPLOT)
  END DO

  PRINT *,' '
  PRINT *,' Successful run with results written to files.'
  PRINT *,' '
  PRINT *,' The data is imported into Matlab as matrices '
  PRINT *,' and S(t) plotted for the 4 runs in one figure by'
  PRINT *,' '
  PRINT *,' >> [t1,y1,t2,y2,t3,y3,t4,y4] = infection;'
  PRINT *,' '
  PRINT *,' The data is imported into Matlab as matrices '
  PRINT *,' and I(t) plotted for the 4 runs in one figure by'
  PRINT *,' '
  PRINT *,' >> [tI1,I1,tI2,I2,tI3,I3,tI4,I4] = infectionIoft;'
  PRINT *,' '

  STOP
END PROGRAM infection
