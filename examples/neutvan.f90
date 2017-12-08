MODULE define_DDEs

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: NEQN=1,NLAGS=1

CONTAINS

  SUBROUTINE DDES(T,Y,Z,DY)
   
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y,DY
!    DOUBLE PRECISION, DIMENSION(NEQN,2*NLAGS) :: Z
    DOUBLE PRECISION, DIMENSION(:,:) :: Z
    DOUBLE PRECISION, DIMENSION(NEQN,NLAGS) :: YLAG,DYLAG
    INTENT(IN)  :: T,Y,Z
    INTENT(OUT) :: DY
    ! NOTE:
    ! This is a neutral problem. Columns 1...NLAGS of the Z
    ! matrix contain the delayed solution values and columns
    ! NLAGS+1...2*NLAGS contain the delayed derivative values.
     
      YLAG  = Z(:,1:NLAGS)
      DYLAG = Z(:,NLAGS+1:2*NLAGS)

      DY(1) = COS(T)*(1D0 + YLAG(1,1)) + 0.3D0*Y(1)*DYLAG(1,1)  &
        + 0.7D0*SIN(T)*COS( T*SIN(T)**2 ) - SIN( T+T*SIN(T)**2 )
                   
    RETURN
  END SUBROUTINE DDES

  SUBROUTINE HISTORY(T,Y)

    DOUBLE PRECISION :: T
    !DOUBLE PRECISION, DIMENSION(NEQN,2) :: Y
    DOUBLE PRECISION, DIMENSION(:) :: Y
    INTENT(IN)  :: T
    INTENT(OUT) :: Y
    ! NOTE:
    ! This is a neutral problem.  Return the history
    ! in column 1 of the Y matrix and the derivative
    ! of the history in column 2.
    ! NEW: 06/28/2013:
    ! Return solution in 1,...,NEQN
    ! Derivative = NEQN+1,...,2*NEQN
    ! Y(1,1) = SIN(T)
    ! Y(1,2) = COS(T)
      Y(1) = SIN(T)
      Y(2) = COS(T)

    RETURN
  END SUBROUTINE HISTORY

  SUBROUTINE BETA(T,Y,BVAL)
    
    DOUBLE PRECISION :: T
    DOUBLE PRECISION, DIMENSION(:) :: Y
    DOUBLE PRECISION, DIMENSION(:) :: BVAL
    INTENT(IN)  :: T,Y
    INTENT(OUT) :: BVAL

      BVAL(1) = T*Y(1)**2

    RETURN
  END SUBROUTINE BETA

  SUBROUTINE CHECK(T,Y,ERROR)
    
    DOUBLE PRECISION :: T,Y
    DOUBLE PRECISION :: ERROR
    INTENT(IN)  :: T,Y
    INTENT(OUT) :: ERROR

      ERROR = ABS(Y - SIN(T))

    RETURN
  END SUBROUTINE CHECK

END MODULE define_DDEs

!******************************************************************

PROGRAM neutvan 

! Program to solve the last example in the SGT book. This 
! is a neutral problem with a state-dependent vanishing delay.
!
! The DDE is defined in the module define_DDEs.  The problem
! is solved here with DDE_SOLVER and its output written to a
! file.  The auxilary function NEUTVAN.M imports the data into
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

  TYPE(DDE_SOL) :: SOL 
  ! The fields of SOL are expressed in terms of the
  ! number of differential equations, NEQN, and the 
  ! number of output points, NPTS:

  !   SOL%NPTS      -- NPTS,number of output points.
  !
  !   SOL%T(NPTS)   -- values of independent variable, T.  
  !
  !   SOL%Y(NPTS,NEQN) -- values of dependent variable, Y,
  !                       corresponding to values of SOL%T. 

  TYPE(DDE_OPTS) :: OPTS

  ! Local variables:
  INTEGER :: I,J
  !DOUBLE PRECISION :: ERROR,MAXERROR
  DOUBLE PRECISION :: ERROR,MAXERROR,TCHECK,YCHECK

  OPTS = DDE_SET(RE=1D-10,AE=1D-10,NEUTRAL=.TRUE.)

  SOL = DDE_SOLVER(NVAR,DDES,BETA,HISTORY,&
                  (/0D0,ASIN(1D0)/),OPTIONS=OPTS)

  ! Was the solver successful?
  IF(SOL%FLAG == 0) THEN

    PRINT *,' '    
    ! Print integration statistics:
    CALL PRINT_STATS(SOL)
 
    ! Write the solution to a file for subsequent 
    ! plotting in Matlab.
    OPEN(UNIT=6, FILE='neutvan.dat')
    DO I = 1,SOL%NPTS
      WRITE(UNIT=6,FMT='(2D12.4)') SOL%T(I),(SOL%Y(I,J),J=1,NEQN)
    ENDDO

    PRINT *,' Normal return from DDE_SOLVER with results'
    PRINT *," written to the file 'neutvan.dat'."
    PRINT *,' '
    PRINT *,' These results can be accessed in Matlab'
    PRINT *,' and plotted by'
    PRINT *,' '
    PRINT *," >> [t,y] = neutvan;"
    PRINT *,' '
    
    MAXERROR = 0D0
    DO I = 1,SOL%NPTS
       ! CALL CHECK(SOL%T(I),SOL%Y(I,1),ERROR)
       TCHECK = SOL%T(I)
       YCHECK = SOL%Y(I,1)
       CALL CHECK(TCHECK,YCHECK,ERROR)      
       MAXERROR = MAX(ERROR,MAXERROR)
    END DO
    PRINT *,' '
    PRINT *,' Comparison to the analytical solution shows'
    PRINT *,' that at all steps the absolute error is less'
    PRINT *,' than ',MAXERROR,'.'
    PRINT *,' '

  ELSE

    PRINT *,' Abnormal return from DDE_SOLVER with FLAG = ',&
    SOL%FLAG

  ENDIF
 
  CALL RELEASE_ARRAYS(SOL,OPTS)
   
  STOP
END PROGRAM neutvan

