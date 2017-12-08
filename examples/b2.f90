MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
   
    DOUBLE PRECISION :: T,YLAG,F
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    INTENT (IN) :: T, Y, Z
    INTENT (OUT) :: DY
     
      ylag = z(1,1)
      if (ylag < 0d0) then
         f = 1d0
      else
         f = -1d0
      endif
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

PROGRAM b2 

! Program to solve ddesd example B2.
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

  INTEGER, DIMENSION(2) :: NVAR = (/NEQN,NLAGS/)

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

  ! Local variables:
  INTEGER :: I,J
  DOUBLE PRECISION :: ERROR,MAXERROR,TCHECK
  DOUBLE PRECISION, DIMENSION(NEQN) :: YCHECK
  CHARACTER(7+6*NEQN) :: EXPORT
  ! FNAME is the name of the output file.
  CHARACTER(10) :: FNAME='export.dat'

  OPTS = DDE_SET(RE=1D-5,AE=1D-5)
  SOL = DDE_SOLVER(NVAR, DDES, BETA, HISTORY, &
                   (/ 0D0, 2D0*LOG(66D0) /), OPTIONS=OPTS)

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

  ELSE

    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG

  ENDIF

  ! Integration statistics:
  CALL PRINT_STATS(SOL)
  
  CALL RELEASE_ARRAYS(SOL,OPTS)

END PROGRAM b2

