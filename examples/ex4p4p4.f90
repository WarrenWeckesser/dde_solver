MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1,NEF=2

! Physical parameters:
  DOUBLE PRECISION, PARAMETER :: R=3.5D0,M=19D0

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN) :: T,Y,Z
    INTENT(OUT) :: DY
   
    DY(1) = R*Y(1)*(1D0 - Z(1,1)/M)

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y
    INTENT(IN) :: T
    INTENT(OUT) :: Y

    ! Use the JUMPS option to tell solver about
    ! the discontinuity at t = 0.
    IF (T == 0D0) THEN
       Y(1) = 19.001D0
    ELSE
       Y(1) = 19D0
    END IF

    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE EF(T,Y,DY,Z,G)

     DOUBLE PRECISION :: T
     DOUBLE PRECISION, DIMENSION(:) :: Y,DY
     DOUBLE PRECISION, DIMENSION(:,:) :: Z
     DOUBLE PRECISION, DIMENSION(:) :: G
     INTENT(IN) :: T,Y,DY,Z
     INTENT(OUT) :: G

     ! Locate extrema as points where the derivative
     ! vanishes.  The DIRECTION option is used to
     ! distinguish maxima and minima.
     G = (/ DY(1),DY(1) /)

     RETURN
  END SUBROUTINE EF


END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p4p4 

! Example 4.4.4 in the SGT book.
!
! The DDE is defined in the module define_DDEs.  The problem
! is solved here with DDE_SOLVER and its output written to a
! file.  The auxilary function EX4P4P4.M imports the data into
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
  TYPE(DDE_INT) :: YINT
  TYPE(DDE_OPTS) :: OPTS

  DOUBLE PRECISION, PARAMETER :: T0=0D0,TFINAL=40D0
  ! Prepare output points
  INTEGER, PARAMETER :: NOUT=1000
  DOUBLE PRECISION, DIMENSION(NOUT) :: TINT

  ! Local variables:
  INTEGER :: I

  OPTS = DDE_SET(RE=1D-4,AE=1D-7,INTERPOLATION=.TRUE.,&
                 JUMPS=(/ 0D0 /),DIRECTION=(/ +1,-1 /))
 
  SOL = DDE_SOLVER(NVAR,DDES,(/ 0.74D0 /),HISTORY, &
         (/ T0,TFINAL /),OPTIONS=OPTS,EVENT_FCN=EF) 

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN

     ! Form values for phase plane plot:
     TINT = (/ (T0+(I-1)*((TFINAL-T0)/(NOUT-1)), I=1,NOUT) /)
     YINT = DDE_VAL(TINT,SOL,DERIVATIVES=.TRUE.)

     ! Write the solution to a file for subsequent plotting
     ! in Matlab.
     OPEN(UNIT=6, FILE='ex4p4p4.dat')
     DO I = 1,NOUT
       WRITE(UNIT=6,FMT='(3D12.4)') TINT(I),YINT%YT(I,1),&
                                    YINT%DT(I,1)
     END DO

     ! Write the extrema to a file for plotting.
     OPEN(UNIT=7,FILE='ex4p4p4extr.dat')
     DO I = 1,SOL%NE
       WRITE(UNIT=7,FMT='(I10,2D12.4)') SOL%IE(I), &
                                        SOL%TE(I),SOL%YE(I,1)
     END DO
    
     PRINT *,' Normal return from DDE_SOLVER with results'
     PRINT *," written to the file 'ex4p4p4.dat' and the"
     PRINT *," extrema written to 'ex4p4p4extr.dat'."
     PRINT *,' '
     PRINT *,' These results can be accessed in Matlab'
     PRINT *,' and plotted by'
     PRINT *,' '
     PRINT *," >> [t,y,yp,te,ye,ie] = ex4p4p4;"
     PRINT *,' '

  ELSE
 
     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG
 
  ENDIF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
   
  CALL RELEASE_ARRAYS(SOL,OPTS)
  CALL RELEASE_INT(YINT)
   
  STOP
END PROGRAM ex4p4p4
