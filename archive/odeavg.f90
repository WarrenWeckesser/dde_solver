  FUNCTION ODEAVG(AVG_OPT,DELAYS,NVAR,ODES,HISTORY,TSPAN,OPTIONS, &
       EVENT_FCN,CHANGE_FCN,OUT_FCN,USER_TRIM_GET) RESULT (SOL)

    ! To use this interface:

    !     AVG_OPT is an array of two entries, [NAVG,KIND].  If NAVG = 1, y(t)
    !     is averaged over an interval of length DELTA=DELAYS(1)>0; and if
    !     NAVG = 2, a second average is computed. If KIND = 1, the moving
    !     average is computed and if KIND = 2, the square of the RMS norm
    !     of the solution is computed.
    !     Summary: AVG_OPT =
    !       <1,1> means the solution will be averaged one time
    !             (once averaged solution)
    !       <2,1> means the solution will be averaged two times
    !             (twice averaged solution)
    !       <1,2> means y^2(t) will be averaged one time
    !             (RMS amplitude)
    !       <2,2> means y^2(t) will be averaged two times
    !             (once averaged RMS amplitude)

    !     DELAYS = vector of length 1 containing the moving average
    !              increment DELTA

    !     ODES(T,Y,DY) = name of subroutine in which the ODE derivatives
    !                    will be defined

    !     HISTORY = vector containing the initial conditions

    !     The remaining arguments have the same meaning as
    !     in the DKL_* drivers.

    !     *****Time shifted solutions are returned by ODEAVG*****
    !     The solution that is returned for, say, TVAL was computed
    !     at the time TVAL + DELINC where DELINC = DELTA/2 if one
    !     average was performed; and DELINC = DELTA if two averages
    !     were performed. DELINT is returned in the solution structure
    !     field SOL%TSHIFT.

    !     Example:
    !     TSPAN(1) =    0 = initial time
    !     TSPAN(2) = 1000 = final time
    !     DELAVG = 25 = moving average increment
    !     AVG_OPT = <1,1>
    !     ODEAVG will return
    !        SOL%T = TSPAN(1)+DELTA/2,...,TSPAN(2)-DELTA/2 = 12.5,...,987.5.
    !     The corresponding once averaged solution values will have
    !     been computed at the times TSPAN(1)+DELTA,...,TSPAN(2)=0,...1000.
    !     Example:
    !     TSPAN(1) =    0 = initial time
    !     TSPAN(2) = 1000 = final time
    !     DELAVG = 25 = moving average increment
    !     AVG_OPT = <2,1>
    !     ODEAVG will return
    !        SOL%T = TSPAN(1),...,TSPAN(2)-DELAVG = 0,...,975.
    !     The corresponding twice averaged solution values will have
    !     been computed at the times DTSPAN(1)+DELAVG,...,TSPAN(2)=25,...,1000.

    !     If you use the INTERPOLATION option and subsequently call
    !     DDE_VAL to obtain the solution for times in an array TINT(*),
    !     pass the values TINT(*)+SOL%TSHIFT to DDE_VAL to obtain
    !     the solution corresponding to TINT(*).

    ! .. Function Return Value ..
    TYPE (DDE_SOL), TARGET :: SOL
    ! ..
    ! .. Structure Arguments ..
    TYPE (DDE_OPTS), OPTIONAL :: OPTIONS
    ! ..
    ! .. Array Arguments ..
    DOUBLE PRECISION, TARGET :: DELAYS(:), HISTORY(:)
    DOUBLE PRECISION :: TSPAN(:)
    INTEGER :: NVAR(:), AVG_OPT(:)
    ! ..
    ! .. Subroutine Arguments ..
    OPTIONAL :: OUT_FCN, CHANGE_FCN, EVENT_FCN, USER_TRIM_GET
    ! ..
    INTERFACE
       SUBROUTINE ODES(T,Y,DY)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
         INTENT(IN):: T,Y
         INTENT(OUT) :: DY
       END SUBROUTINE ODES
    END INTERFACE
    INTERFACE
       SUBROUTINE CHANGE_FCN(NEVENT,TEVENT,YEVENT,DYEVENT,HINIT, &
            DIRECTION,ISTERMINAL,QUIT)
         INTEGER :: NEVENT
         INTEGER, DIMENSION(:) :: DIRECTION
         DOUBLE PRECISION :: TEVENT,HINIT
         DOUBLE PRECISION, DIMENSION(:) :: YEVENT,DYEVENT
         LOGICAL :: QUIT
         LOGICAL, DIMENSION(:) :: ISTERMINAL
         INTENT(IN) :: NEVENT,TEVENT
         INTENT(INOUT) :: YEVENT,DYEVENT,HINIT,DIRECTION,ISTERMINAL,QUIT
       END SUBROUTINE CHANGE_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE EVENT_FCN(T,Y,DYDT,Z,G)
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DYDT
         DOUBLE PRECISION, DIMENSION(:,:) :: Z
         DOUBLE PRECISION, DIMENSION(:) :: G
         INTENT(IN):: T,Y,DYDT,Z
         INTENT(OUT) :: G
       END SUBROUTINE EVENT_FCN
    END INTERFACE
    INTERFACE
       SUBROUTINE OUT_FCN(T,Y,DY,N,NEVENT)
         INTEGER :: N,NEVENT
         DOUBLE PRECISION :: T
         DOUBLE PRECISION, DIMENSION(:) :: Y,DY
       END SUBROUTINE OUT_FCN
    END INTERFACE
   INTERFACE
      SUBROUTINE USER_TRIM_GET
      END SUBROUTINE USER_TRIM_GET
   END INTERFACE
    ! ..
    ! .. Local Structures ..
    TYPE (DDE_OPTS), TARGET :: OPTS
    ! ..
    ! .. Local Scalars ..
    INTEGER :: I, IER, IFLAG, NEQN, NGUSER, NLAGS, NOUT, NPTS, &
         KFIRST, NPTS_TEMP, MKIND, IQUEUE, IDROPM1, ISAVE
    DOUBLE PRECISION :: TF_TEMP, T0_TEMP
    DOUBLE PRECISION, ALLOCATABLE :: TEMP(:), YTEMP(:,:)
    ! ..
    ! .. Intrinsic Functions ..
    INTRINSIC PRESENT, SIZE
    ! ..
    ! .. User subroutines ..
    !   EXTERNAL ODES,OUT_FCN,CHANGE_FCN,EVENT_FCN,USER_TRIM_GET
    ! ..
    SOLIPOINT_IS_ALLOCATED = .FALSE.

    SHIFTT = .TRUE.
    SOL%SHIFT = SHIFTT

    HAVE_TRIM_GET = PRESENT(USER_TRIM_GET)

    IF (PRESENT(OPTIONS)) THEN
       OPTS = OPTIONS
    ELSE
       OPTS = DDE_SET()
    END IF

    MAVG = AVG_OPT(1)
    MKIND = AVG_OPT(2)
    IF (MAVG<1 .OR. MAVG>2) THEN
       PRINT *, ' You must input the number of averages, AVG_OPT(1).'
       PRINT *, ' Legal values are 1 and 2.'
       STOP
    END IF
    IF (MKIND<1 .OR. MKIND>2) THEN
       PRINT *, ' You must input the kind of average, AVG_OPT(2).'
       PRINT *, ' Legal values are 1 and 2.'
       STOP
    END IF
    AVG_TYPE = MKIND
    IF (MAVG==1) THEN
       DELINC = 0.5D0 * DELAVG
    ELSE
       DELINC = DELAVG
    END IF

    IF (SIZE(NVAR)<1) THEN
       PRINT *, ' You must input the number of equations, NVAR(1).'
       STOP
    ELSE
       DELAVG = DELAYS(1)
       DROPZ = .TRUE.
       NEQN = NVAR(1)
       NLAGS = 1
       NEQN_USER = NEQN
       NLAGS_USER = NLAGS
    END IF

    SOL%TSHIFT = DELINC
    NEQN = (MAVG + 1) * NEQN
    HAVE_EVENT_FCN = PRESENT(EVENT_FCN)
    IF (HAVE_EVENT_FCN) THEN
       IF (SIZE(NVAR)<3) THEN
          PRINT *, ' You must input the number of event functions, NVAR(3).'
          STOP
       ELSE
          NGUSER = NVAR(3)
       END IF
       ! Provide for 10 events.
       ALLOCATE (SOL%IE(10),SOL%TE(10),SOL%YE(10,NEQN),STAT=IER)
       MYIPOINT(15) = 1
       MYIPOINT(10) = 1
       MYIPOINT(11) = 1
       CALL CHECK_STAT(IER,13)
    ELSE
       NGUSER = 0
    END IF
    SOL%NE = 0
    CALL EXPAND_OPTS(NEQN,NGUSER,OPTS)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1204)
    ! END IF
    ! IF (OPTS%NEUTRAL) THEN
    !    ALLOCATE(YOLD(2*NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1205)
    ! ELSE
    !    ALLOCATE(YOLD(NEQN),STAT=IER)
    !    CALL CHECK_STAT(IER,1206)
    !    PRINT *, ' YOLD has been allocated.'
    ! END IF

    ! The following was moved to DDE_DRV1.
    ! YOLD(1:NVAR(1)) = HISTORY
    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    NSPAN = SIZE(TSPAN)
    IF (HAVE_OUT_FCN .AND. NSPAN>2) THEN
       ! Change to output at every step.
       TSPAN = (/ TSPAN(1), TSPAN(NSPAN) /)
       NSPAN = 2
    END IF
    CALL CHECK_TSPAN(TSPAN,NOUT)
    ALLOCATE (SOL%T(NOUT),SOL%Y(NOUT,NEQN),STAT=IER)
    MYIPOINT(8) = 1
    MYIPOINT(9) = 1
    CALL CHECK_STAT(IER,14)
    SOL%NPTS = 0
    PASS_SOL => SOL
    NULLIFY (PASS_DELAYS,PASS_HISTORY)
    PASS_DELAYS => DELAYS
    CONSTANT_DELAYS = .TRUE.
    PASS_HISTORY => HISTORY
    CONSTANT_HISTORY = .TRUE.
    HAVE_OUT_FCN = PRESENT(OUT_FCN)
    IF (PRESENT(CHANGE_FCN)) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,CHANGE_FCN, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .FALSE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,CHANGE_FCN, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE IF (HAVE_EVENT_FCN) THEN
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,EVENT_FCN,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    ELSE
       IF (HAVE_OUT_FCN) THEN
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,DUMMY_GUSER,DUMMY_CHANGE, &
                  OUT_FCN,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       ELSE
          IF (HAVE_TRIM_GET) THEN
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,USER_TRIM_GET)
          ELSE
             DONT_CALL_CHANGE = .TRUE.
             CALL DDE_DRV1(NEQN,NLAGS,OPTS,ODES,DUMMY_GUSER,DUMMY_CHANGE, &
                  SOL_OUT,DUMMY_BETA,DUMMY_YINIT,TSPAN,NGUSER,TRIM_GET_DUMMY)
          END IF
       END IF
    END IF

    IFLAG = ERROR_FLAG

    ! Generate an error message if the integration was not successful.
    IF (IFLAG/=0) THEN
       PRINT *, ' One or more errors occurred in ODEAVG.'
    END IF

    ! Save the necessary interpolation structure if necessary.
    SOLIPOINT_IS_ALLOCATED = .FALSE.
    IF (IFLAG==0 .AND. OPTS%INTERPOLATION) THEN
       SOL%SHIFT = SHIFTT
       ! Add the interpolation arrays to SOL.
       ! Reset the dimensions of TQUEUE and QUEUE to return only
       ! the portions of TQUEUE and QUEUE which were actually used
       ! (okay to reset since the arrays are about to be deallocated).
       IF (AVG_TYPE==0) THEN
          LTQUEUE = MYIPOINT(2)
          MYQCOLS = 10*LTQUEUE
          MYIPOINT(1) = MYQCOLS
          ALLOCATE (SOL%IPOINT(LIPOINT),SOL%YOFT(NEQN), &
               SOL%QUEUE(NEQN,MYQCOLS), &
               SOL%TQUEUE(0:LTQUEUE),STAT=IER)
          CALL CHECK_STAT(IER,15)
          SOLIPOINT_IS_ALLOCATED = .TRUE.
          MYIPOINT(16) = 1
          MYIPOINT(13) = 1
          MYIPOINT(12) = 1
          MYIPOINT(14) = 1
          SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
          SOL%QUEUE(1:NEQN,1:MYQCOLS) = &
               QUEUE(1:NEQN,1:MYQCOLS)
          SOL%TQUEUE(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
       ELSE
          !    Return only the most averaged solution.
          LTQUEUE = MYIPOINT(2)
          MYQCOLS = 10*LTQUEUE
          MYIPOINT(1) = MYQCOLS
          ALLOCATE (SOL%IPOINT(LIPOINT),SOL%YOFT(NEQN_USER), &
               SOL%QUEUE(NEQN_USER,MYQCOLS), &
               SOL%TQUEUE(0:LTQUEUE),STAT=IER)
          CALL CHECK_STAT(IER,15)
          SOLIPOINT_IS_ALLOCATED = .TRUE.
          MYIPOINT(16) = 1
          MYIPOINT(13) = 1
          MYIPOINT(12) = 1
          MYIPOINT(14) = 1
          SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
          SOL%QUEUE(1:NEQN_USER,1:MYQCOLS) = &
               QUEUE(NEQN-NEQN_USER+1:NEQN,1:MYQCOLS)
          SOL%TQUEUE(0:LTQUEUE) = TQUEUE(0:LTQUEUE)
       END IF
    END IF

    IF (HAVE_TRIM_GET) THEN
       ! Write the remaining queue information.
       IQUEUE = MYIPOINT(2)
       DO I = 1, IQUEUE
          ISAVE = I
          IF (TQUEUE(I)>TQLAST) THEN

             IF (DEBUG) THEN
                WRITE(DBGUNIT,9999)
                WRITE(DBGUNIT,9998)
                9999 FORMAT(' About to call TRIM_SAVE following the')
                9998 FORMAT(' completion of the integration.')
             END IF

             CALL TRIM_SAVE(I,IQUEUE)
             CALL USER_TRIM_GET
             GOTO 10
          END IF
       END DO
10     CONTINUE

       IF (DEBUG) THEN
          IDROPM1 = IQUEUE - ISAVE + 1
          WRITE(DBGUNIT,9997) TQUEUE(IQUEUE)
          WRITE(DBGUNIT,9996) IDROPM1
          9997 FORMAT(' TNEW = ', D20.10)
          9996 FORMAT(' Dropped from end of queue = ', I10)
       END IF

       TQLAST = TQUEUE(IQUEUE)

    END IF

    ! Deallocate the remaining DDE_SOLVER arrays.
    DEALLOCATE (QUEUE,TQUEUE,STAT=IER)
    CALL CHECK_STAT(IER,16)

    ! IF (ALLOCATED(YOLD)) THEN
    !    DEALLOCATE(YOLD,STAT=IER)
    !    CALL CHECK_STAT(IER,1250)
    ! END IF

    ! Derivative sign change array.
    ! IF (AVG_TYPE /= 0) THEN
    !    DO I = 1, NEQN_USER
    !       IF (DYSIGNS(NEQN-NEQN_USER+I) /= 0) THEN
    !          PRINT *, ' Component ', I , ' of the most averaged solution'
    !          PRINT *, ' changed sign ', DYSIGNS(NEQN-NEQN_USER+I), ' times.'
    !       END IF
    !    END DO
    !    DEALLOCATE (DYSIGNS,STAT=IER)
    !    CALL CHECK_STAT(IER,4)
    ! END IF

    ! Trim solution structure for output.
    CALL CONTRACT_SOL(PASS_SOL,NEQN)

    ! Drop the solution before DELAVG and shift the
    ! independent variable.

    IF (SHIFTT .AND. (MAVG==1 .OR. MAVG==2)) THEN
       NPTS = PASS_SOL%NPTS
       T0_TEMP = PASS_SOL%T(1)
       TF_TEMP = PASS_SOL%T(NPTS)
       DO I = 1, NPTS
          IF (MAVG==1 .AND. PASS_SOL%T(I)>PASS_SOL%T(1)+DELAVG) THEN
             KFIRST = I - 1
             GOTO 20
          END IF
          IF (MAVG==2 .AND. PASS_SOL%T(I)>PASS_SOL%T(1)+2.0D0*DELAVG) THEN
             KFIRST = I - 1
             GOTO 20
          END IF
       END DO
       PRINT *, ' An error occurred in ODEAVG during the shift'
       PRINT *, ' of the independent variable.'
       RETURN
20     CONTINUE
       NPTS_TEMP = NPTS - (KFIRST-1)
       IF (MAVG==2) THEN
          DO I = 1, NPTS_TEMP
             PASS_SOL%T(I) = PASS_SOL%T(KFIRST-1+I) - DELAVG
          END DO
       END IF
       IF (MAVG==1) THEN
          DO I = 1, NPTS_TEMP
             PASS_SOL%T(I) = PASS_SOL%T(KFIRST-1+I) - 0.5D0*DELAVG
          END DO
       END IF
       PASS_SOL%T(1) = MAX(T0_TEMP,PASS_SOL%T(1))
       PASS_SOL%T(NPTS_TEMP) = MIN(TF_TEMP,PASS_SOL%T(NPTS_TEMP))
       DO I = 1, NPTS_TEMP
          PASS_SOL%Y(I,1:NEQN_USER) = PASS_SOL%Y(KFIRST-1+I,1:NEQN_USER)
       END DO
       ALLOCATE (TEMP(NPTS_TEMP),STAT=IER)
       CALL CHECK_STAT(IER,38)
       TEMP = PASS_SOL%T(1:NPTS_TEMP)
       DEALLOCATE (PASS_SOL%T,STAT=IER)
       CALL CHECK_STAT(IER,39)
       ALLOCATE (PASS_SOL%T(NPTS_TEMP),STAT=IER)
       MYIPOINT(8) = 1
       CALL CHECK_STAT(IER,40)
       PASS_SOL%T = TEMP
       DEALLOCATE (TEMP,STAT=IER)
       CALL CHECK_STAT(IER,41)
       ALLOCATE (YTEMP(NPTS_TEMP,NEQN_USER),STAT=IER)
       CALL CHECK_STAT(IER,42)
       YTEMP = PASS_SOL%Y(1:NPTS_TEMP,1:NEQN_USER)
       DEALLOCATE (PASS_SOL%Y,STAT=IER)
       CALL CHECK_STAT(IER,43)
       ALLOCATE (PASS_SOL%Y(NPTS_TEMP,NEQN_USER),STAT=IER)
       MYIPOINT(9) = 1
       PASS_SOL%Y = YTEMP
       CALL CHECK_STAT(IER,44)
       DEALLOCATE (YTEMP,STAT=IER)
       CALL CHECK_STAT(IER,44)
       PASS_SOL%NPTS = NPTS_TEMP
    END IF

    ALLOCATE (SOL%STATS(6),STAT=IER)
    CALL CHECK_STAT(IER,3)
    MYIPOINT(17) = 1
    IF (SOLIPOINT_IS_ALLOCATED) THEN
    ELSE
       ALLOCATE (SOL%IPOINT(LIPOINT),STAT=IER)
       CALL CHECK_STAT(IER,3)
       SOLIPOINT_IS_ALLOCATED = .TRUE.
       MYIPOINT(16) = 1
    END IF
    SOL%IPOINT(1:LIPOINT) = MYIPOINT(1:LIPOINT)
    SOL%FLAG = IFLAG
    SOL%STATS(1) = NSTEPS
    SOL%STATS(2) = NFAILS
    SOL%STATS(3) = NFEVAL
    SOL%STATS(4) = NFAILC
    SOL%STATS(5) = ARRAY_STORAGE
    SOL%STATS(6) = ROOT_FUNCTIONS

    MYIPOINT(1:LIPOINT) = 0

    IF (ASSOCIATED(PASS_SOL)) NULLIFY(PASS_SOL)
    IF (ASSOCIATED(PASS_DELAYS)) NULLIFY(PASS_DELAYS)
    IF (ASSOCIATED(PASS_HISTORY)) NULLIFY(PASS_HISTORY)

    RETURN
  END FUNCTION ODEAVG
