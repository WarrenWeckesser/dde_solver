MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=4,NLAGS=1,NEF=1

  INTEGER :: STATE

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)

    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN):: T,Y,Z
    INTENT(OUT) :: DY
    DOUBLE PRECISION, DIMENSION(NEQN) :: YLAG

    ! Physical parameters
    DOUBLE PRECISION, PARAMETER :: H1=2D0,H2=0.8D0,&
    H3=1D4,H4=0.17D0,H5=0.5D0,H6=3D2,H7=0.12D0,H8=8D0
    DOUBLE PRECISION :: XI

      YLAG = Z(:,1)

      IF (STATE==1) THEN
         XI = 1D0
      ELSE
         XI = (10D0/9D0)*(1D0 - Y(4))
      END IF

      DY(1) = (H1 - H2*Y(3))*Y(1)
      DY(2) = XI*H3*YLAG(3)*YLAG(1) - H5*(Y(2) - 1D0)
      DY(3) = H4*(Y(2) - Y(3)) - H8*Y(3)*Y(1)
      DY(4) = H6*Y(1) - H7*Y(4)

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)

     DOUBLE PRECISION :: T
     DOUBLE PRECISION, DIMENSION(:) :: Y
     INTENT(IN):: T
     INTENT(OUT) :: Y

       Y(1) = MAX(0D0,T+1D-6)
       Y(2) = 1D0
       Y(3) = 1D0
       Y(4) = 0D0

    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE EF(T,Y,DY,Z,G)

     DOUBLE PRECISION :: T
     DOUBLE PRECISION, DIMENSION(:) :: Y,DY
     DOUBLE PRECISION, DIMENSION(:,:) :: Z
     DOUBLE PRECISION, DIMENSION(:) :: G
     INTENT(IN):: T,Y,DY,Z
     INTENT(OUT) :: G

       G(1) = Y(4) - 1D-1

     RETURN
  END SUBROUTINE EF

  SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
                  DIRECTION,ISTERMINAL,QUIT)

    INTEGER :: NEVENT
    LOGICAL, DIMENSION(:) :: ISTERMINAL
    INTEGER, DIMENSION(:) :: DIRECTION
    DOUBLE PRECISION :: TEVENT,HINIT
    DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
    ! NOTE:  THE INPUT VALUE OF QUIT IS .FALSE., SO IT NEED BE
    !        SET ONLY IF THE INTEGRATION IS TO BE TERMINATED.
    LOGICAL :: QUIT
    INTENT(IN) :: NEVENT,TEVENT
    INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT

      STATE = -STATE
      HINIT = 1D-7

    RETURN
  END SUBROUTINE CHNG

END MODULE define_DDEs

!******************************************************************

PROGRAM ex4p8 

  ! Marchuk immunology model, Exercise 4.8 in the book.

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

  DOUBLE PRECISION, DIMENSION(NLAGS) :: LAGS=(/ 0.5D0 /)
  DOUBLE PRECISION, DIMENSION(2) :: TSPAN=(/ 0D0,60D0 /)

  TYPE(DDE_OPTS) :: OPTS

  ! Local variables:
  INTEGER :: I,J

  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(10) :: FNAME='export.dat'  

  STATE = 1

  OPTS = DDE_SET(RE=1D-5,AE=1D-8,JUMPS=(/ -1D-6 /))

  SOL = DDE_SOLVER(NVAR,DDES,LAGS,HISTORY,TSPAN,&
           OPTIONS=OPTS,EVENT_FCN=EF,CHANGE_FCN=CHNG) 

  ! Was the solver successful?
  IF (SOL%FLAG == 0) THEN
    
     ! Print integration statistics:
     CALL PRINT_STATS(SOL)
 
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
     PRINT *,' '
     PRINT *,' To access the output in Matlab and plot the '
     PRINT *,' solution components scaled as in Fig. 15.8 '
     PRINT *,' of the Hairer, Norsett, and Wanner book, use'
     PRINT *,' '
     PRINT *,' >> [t,y] = ex4p8;'
     PRINT *,' '

  ELSE

     PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
     SOL%FLAG

  ENDIF
  
  STOP
END PROGRAM ex4p8
