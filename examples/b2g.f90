MODULE define_DDEs

! Modified version of B2.

  IMPLICIT NONE
  
  INTEGER STATE
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1,NEF=1

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
   
    DOUBLE PRECISION :: T,F
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN) :: T, Y, Z
    INTENT(OUT) :: DY
     
      if (state == 1) then
         f = -1d0
      else
         f = 1d0
      end if
      dy(1) = f - y(1)

    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)

    DOUBLE PRECISION, INTENT(IN) :: T
    DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: Y

      Y(1) = 1D0

    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE BETA(T,Y,BVAL)
    
    DOUBLE PRECISION, INTENT(IN) :: T
    DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: Y
    DOUBLE PRECISION, DIMENSION(:), INTENT(OUT) :: BVAL

      BVAL(1) = T / 2D0

    RETURN
  END SUBROUTINE BETA
  SUBROUTINE EF(T,Y,DY,Z,G)

    DOUBLE PRECISION :: T,YLAG
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:) :: G
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT(IN) :: T, Y, DY, Z
    INTENT(OUT) :: G

      YLAG = Z(1,1)
      G(1) = YLAG
    
    RETURN
  END SUBROUTINE EF

 SUBROUTINE CHNG(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
     DIRECTION,ISTERMINAL,QUIT)

     INTEGER :: NEVENT
     LOGICAL, DIMENSION(:) :: ISTERMINAL
     INTEGER, DIMENSION(:) :: DIRECTION
     DOUBLE PRECISION :: TEVENT,HINIT
     DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
     LOGICAL :: QUIT
     INTENT(IN) :: NEVENT, TEVENT
     INTENT(INOUT) :: YEVENT, DYEVENT, HINIT, DIRECTION, ISTERMINAL, QUIT

       IF (NEVENT == 1) THEN
          ! Reset the derivative switch for ylag_1(t).
          STATE = -STATE
       ENDIF 

    RETURN
  END SUBROUTINE CHNG

  SUBROUTINE CHECK(T,Y,ERROR)
    
    DOUBLE PRECISION :: T,V
    DOUBLE PRECISION, DIMENSION(NEQN) :: Y
    DOUBLE PRECISION :: ERROR

      if (t <= 2d0*log(2d0)) then
          v = 2d0*exp(-t) - 1d0
      else
          if (t <= 2d0*log(6d0)) then
             v = 1d0 - 6d0*exp(-t);
          else
             if (t <= 2d0*log(66d0)) then
                v = 66d0*exp(-t) - 1d0;
              else
                 print *,'Analytical solution not implemented for t>2*log(66).'
                 stop
             end if
          end if
       end if
       error = abs(y(1)-v)

    RETURN
  END SUBROUTINE CHECK

END MODULE define_DDEs

!******************************************************************

PROGRAM b2g

! Program to solve a modified version of HENS B2. In this version
! root finding is used to locate the points at which ylag_(t) = 0.
! A derivative switch is used to minimize the effect of
! discontinuities in y'(t).
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

  !   SOL%NPTS      -- NPTS,number of output points.
  !
  !   SOL%T(NPTS)   -- values of independent variable, T.  
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                 corresponding to values of SOL%T. 

  TYPE(DDE_OPTS) :: OPTS

  LOGICAL, DIMENSION(NEF) :: TERMINAL= (/ .FALSE. /)
  INTEGER, DIMENSION(NEF) :: CROSS_AXIS=(/ 0 /)

  ! Local variables:
  INTEGER :: I,J
  DOUBLE PRECISION :: ERROR,MAXERROR,TCHECK
  DOUBLE PRECISION, DIMENSION(NEQN) :: YCHECK

  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(10) :: FNAME='b2g.dat'  

  ! Format for writing events to an output file.
  CHARACTER(10+6*NEQN) :: EOUT
  ! EFNAME is the name of the output file.
  CHARACTER(13) :: EFNAME='b2gevents.dat'

! Initialize the derivative switch for ylag(t):
  STATE = 1

  OPTS = DDE_SET(ISTERMINAL=TERMINAL,DIRECTION=CROSS_AXIS,&
         RE=1D-5,AE=1D-5)

  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,(/0D0,2D0*LOG(66D0)/),&
        EVENT_FCN=EF,OPTIONS=OPTS,CHANGE_FCN=CHNG)

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
        PRINT *,' >> [t,y,te,ye,ie] = eventsex2D;'
        PRINT *,' '
     END IF

!   Check the accuracy of the computed solution.
    MAXERROR = 0D0
    DO I = 1,SOL%NPTS
       ! CALL CHECK(SOL%T(I),SOL%Y(I,1),ERROR)
       TCHECK = SOL%T(I)
       DO J = 1,NEQN
          YCHECK(J) = SOL%Y(I,J)
       END DO
       CALL CHECK(TCHECK,YCHECK,ERROR)      
       MAXERROR = MAX(ERROR,MAXERROR)
    END DO
    PRINT *,' '
    PRINT *,' Comparison to the analytical solution shows that'
    PRINT *,' at all steps the absolute error is less than ',MAXERROR
    PRINT *,' '

!   Check if the correct number of roots was located.
    IF (SOL%NE > 3 .OR. SOL%NE < 2) THEN
       PRINT *, ' The correct number of roots was not located.'
    ENDIF

!   Check the accuracy of the computed roots.
!   ylag(t) vanishes at 2*log(2), 2*log(6), and 2*log(66).
!   Note: If the numerical solution is negative near TFINAL=2*log(66),
!   a root may be found since, in fact, ylag_1 = 0 there, even though
!   the analytical solution does not change sign.
        DO I = 1,SOL%NE
            IF  (I == 1) ERROR = SOL%TE(1) - 2d0*log(2d0) 
            IF  (I == 2) ERROR = SOL%TE(2) - 2d0*log(6d0) 
            IF  (I == 3) ERROR = SOL%TE(3) - 2d0*log(66d0) 
            PRINT *, ' Error in computed root = ', ERROR
        END DO

  ELSE

    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG

  ENDIF

  CALL RELEASE_ARRAYS(SOL,OPTS)
   
  STOP

END PROGRAM b2g

