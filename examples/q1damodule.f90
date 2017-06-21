MODULE Q1DAMODULE

  IMPLICIT NONE
  
CONTAINS

      SUBROUTINE Q1DA(F, A, B, EPS, R, E, KF, IFLAG, RPAR, IPAR)
!     MODIFIED VERSION OF SUBROUTINE Q1DA.
!       Q1DA IS A SUBROUTINE FOR THE AUTOMATIC EVALUATION
!       OF THE DEFINITE INTEGRAL OF A USER DEFINED FUNCTION
!       OF ONE VARIABLE.
!       From the book "Numerical Methods and Software"
!       by  D. Kahaner, C. Moler, S. Nash
!       Prentice Hall 1988
!       ARGUMENTS IN THE CALL SEQUENCE:
!       A,B   (INPUT) THE ENDPOINTS OF THE INTEGRATION INTERVAL
!       EPS   (INPUT) THE ACCURACY TO WHICH YOU WANT THE INTEGRAL
!             COMPUTED.  IF YOU WANT 2 DIGITS OF ACCURACY SET
!             EPS=.01, FOR 3 DIGITS SET EPS=.001, ETC.
!             EPS MUST BE POSITIVE.
!       R     (OUTPUT) Q1DA'S BEST ESTIMATE OF YOUR INTEGRAL
!       E     (OUTPUT) AN ESTIMATE OF DABS(INTEGRAL-R)
!       KF    (OUTPUT) THE COST OF THE INTEGRATION, MEASURED IN
!             NUMBER OF EVALUATIONS OF YOUR INTEGRAND.
!             KF WILL ALWAYS BE AT LEAST 30.
!       IFLAG (OUTPUT) TERMINATION FLAG...POSSIBLE VALUES ARE
!               0   NORMAL COMPLETION, E SATISFIES
!                        E<EPS  AND  E<EPS*DABS(R)
!               1   NORMAL COMPLETION, E SATISFIES
!                        E<EPS, BUT E>EPS*ABS(R)
!               2   NORMAL COMPLETION, E SATISFIES
!                        E<EPS*ABS(R), BUT E>EPS
!               3   NORMAL COMPLETION BUT EPS WAS TOO SMALL TO
!                     SATISFY ABSOLUTE OR RELATIVE ERROR REQUEST.
!               4   ABORTED CALCULATION BECAUSE OF SERIOUS ROUNDING
!                     ERROR.  PROBABLY E AND R ARE CONSISTENT.
!               5   ABORTED CALCULATION BECAUSE OF INSUFFICIENT STORAGE.
!                     R AND E ARE CONSISTENT.
!               6   ABORTED CALCULATION BECAUSE OF SERIOUS DIFFICULTIES
!                     MEETING YOUR ERROR REQUEST.
!               7   ABORTED CALCULATION BECAUSE EPS WAS SET <= 0.0
!            NOTE: IF IFLAG=3, 4, 5 OR 6 CONSIDER USING Q1DAX INSTEAD.
!    WHERE IS YOUR INTEGRAND?
!    YOU MUST WRITE A FORTRAN FUNCTION, CALLED F, TO EVALUATE
!    THE INTEGRAND.  USUALLY THIS LOOKS LIKE...
!    FUNCTION F(X)
!       F=(EVALUATE THE INTEGRAND AT THE POINT X)
!       RETURN
!    END
!
      IMPLICIT NONE
      DOUBLE PRECISION A, B, E, EPS, FMAX, FMIN, R, W(1000,6), &
      RPAR
      INTEGER KF,IFLAG,NMAX,NINT,IPAR
      DIMENSION RPAR(1),IPAR(1)
      LOGICAL RST
      EXTERNAL F
!
      NINT = 1
      RST = .FALSE.
      NMAX = 1000
      CALL Q1DAX (F, A, B, EPS, R, E, NINT, RST, W, NMAX, FMIN, FMAX,   &
      KF, IFLAG, RPAR, IPAR)
      RETURN
      END SUBROUTINE Q1DA
      SUBROUTINE Q1DAX (F, A, B, EPS, R, E, NINT, RST, W, NMAX, FMIN,   &
      FMAX, KF, IFLAG, RPAR, IPAR)
!
      IMPLICIT NONE
      INTEGER C,NMAX,MXTRY,NINT,KF,IFLAG,I,IROFF,LOC,IPAR
      DOUBLE PRECISION A, B, E, EB, EPMACH, EPS, FMIN, FMINL, FMINR,    &
      FMN, FMX, R, RAB, RABS, RAV, T, TE, TE1, TE2, TR, TR1, TR2, UFLOW,&
      W(NMAX,6), XM, FMAX, FMAXL, FMAXR, RPAR
      DIMENSION RPAR(1),IPAR(1)
      EXTERNAL F
      LOGICAL RST

      EPMACH = 2.22D-16
      UFLOW = 2.23D-308

      MXTRY = NMAX / 2
!     In case there is no more room, we can toss out easy intervals
!     at most MXTRY times.
!     IF (A.EQ.B) THEN
      IF (DABS(B-A).LE.0D0) THEN
         R = 0.D0
         E = 0.D0
         NINT = 0
         IFLAG = 0
         KF = 1
         CALL F(A,RPAR,IPAR,FMIN)
         FMAX = FMIN
         GOTO 20
      ENDIF
      IF (RST) THEN
         IF (IFLAG.LT.3) THEN
            EB = MAX(100.D0*UFLOW,MAX(EPS,50.D0*EPMACH) * DABS(&
            R))
            DO 19 I = 1, NINT
               IF (DABS(W(I,3)) .GT. (EB * (W(I,2) - W(I,1))   &
               / (B - A) ) ) THEN
                  W (I, 3) = DABS(W(I,3) )
               ELSE
                  W (I, 3) = - DABS(W(I,3) )
               ENDIF
   19       END DO
            GOTO 15
         ELSE
            GOTO 20
         ENDIF
      ENDIF
      KF = 0
      IF (EPS.LE.0.D0.OR.NINT.LE.0.OR.NINT.GE.NMAX) THEN
         IFLAG = 7
         GOTO 20
      ENDIF
      IF (NINT.EQ.1) THEN
         W (1,1) = A
         W (2,2) = B
         W (1,5) = A
         W (1,6) = B
         W (2,5) = A
         W (2,6) = B
         W (1,2) = A + (B - A) / 2.D0
         W (2,1) = W (1,2)
         NINT = 2
      ELSE
!        IF (W(1,1).NE.A .OR. W(NINT+1,1).NE.B) THEN
         IF (DABS(W(1,1)-A).NE.0D0 .OR. DABS(W(NINT+1,1)-B).NE.0D0) THEN
            IFLAG = 8
            GOTO 20
         ENDIF
         W (1,5) = A
         DO 89 I = 1, NINT
            W (I,2) = W (I+1,1)
            W (I,5) = W (I,1)
            W (I,6) = W (I,2)
   89    END DO
      ENDIF
!
      IFLAG = 0
      IROFF = 0
      RABS = 0.0D0
      DO 3 I = 1, NINT
         CALL GL15T (F, W(I,1), W(I,2), DBLE(W(I,5)), DBLE( &
         W(I,6)), W(I,4), W(I,3), RAB, RAV, FMN, FMX,       &
         RPAR, IPAR)
         KF = KF + 15
         IF (I.EQ.1) THEN
            R = W(I,4)
            E = W(I,3)
            RABS = RABS + RAB
            FMIN = FMN
            FMAX = FMX
         ELSE
            R = R + W(I,4)
            E = E+W(I,3)
            RABS = RABS + RAB
            FMAX = MAX(FMAX,FMX)
            FMIN = MIN(FMIN,FMN)
         ENDIF
    3 END DO
      DO 10 I = NINT + 1, NMAX
         W (I,3) = 0.D0
   10 END DO
   15 CONTINUE
!
!   MAIN SUBPROGRAM LOOP
!
      IF (100.D0 * EPMACH * RABS.GE.DABS(R).AND.E.LT.EPS) GOTO 20
      EB = MAX(100.D0*UFLOW,MAX(EPS,50.D0*EPMACH)*DABS(R))
      IF (E.LE.EB) GOTO 20
      IF (NINT.LT.NMAX) THEN
         NINT = NINT + 1
         C = NINT
      ELSE
         C = 0
   16    IF (C.EQ.NMAX.OR.MXTRY.LE.0) THEN
            IFLAG = 5
            GOTO 20
         ENDIF
         C = C + 1
         IF (W(C,3).GT.0.0D0) GOTO 16
!            Found an interval to throw out
         MXTRY = MXTRY - 1
      ENDIF
      CALL DAMAX(NINT,W(1,3),1,LOC)
      XM = W (LOC,1) + (W(LOC,2) - W (LOC,1)) / 2.D0
      IF ((MAX(DABS(W(LOC,1)) , DABS(W(LOC,2)))) .GT. ((1. &
      + 100.D0 * EPMACH) * (DABS(XM) + 0.1D+04 * UFLOW))) THEN
         CALL GL15T(F, W(LOC,1), XM, DBLE (W(LOC,5)), DBLE(W(&
         LOC,6)), TR1, TE1, RAB, RAV, FMINL, FMAXL,          &
         RPAR,IPAR)
         KF = KF + 15
         IF (TE1.LT. (EB * (XM - W (LOC, 1) ) / (B - A) ) ) TE1 =      &
         - TE1
         CALL GL15T (F, XM, W(LOC,2), DBLE(W(LOC,5)), DBLE(W(&
         LOC,6)) , TR2, TE2, RAB, RAV, FMINR, FMAXR,                   &
         RPAR,IPAR)
         KF = KF + 15
         FMIN = MIN(FMIN,FMINL,FMINR)
         FMAX = MAX(FMAX,FMAXL,FMAXR)
         IF (TE2.LT. (EB * (W(LOC,2) - XM) / (B - A) ) ) TE2 =        &
         - TE2
         TE = DABS (W(LOC,3))
         TR = W (LOC,4)
         W (C, 3) = TE2
         W (C, 4) = TR2
         W (C, 1) = XM
         W (C, 2) = W(LOC,2)
         W (C, 5) = W(LOC,5)
         W (C, 6) = W(LOC,6)
         W (LOC,3) = TE1
         W (LOC,4) = TR1
         W (LOC,2) = XM
         E = E-TE+ (DABS(TE1) + DABS(TE2) )
         R = R - TR + (TR1 + TR2)
         IF (DABS(DABS(TE1) + DABS(TE2) - TE) .LT.0.001D0 * TE) THEN
            IROFF = IROFF + 1
            IF (IROFF.GE.10) THEN
               IFLAG = 4
               GOTO 20
            ENDIF
         ENDIF
      ELSE
         IF (EB.GT.W(LOC,3)) THEN
            W(LOC,3) = 0.D0
         ELSE
            IFLAG = 6
            GOTO 20
         ENDIF
      ENDIF
      GOTO 15
!     ALL EXITS FROM HERE
   20 CONTINUE
      IF (IFLAG.GE.4) RETURN
      IFLAG = 3
      T = EPS * DABS (R)
      IF (E.GT.EPS.AND.E.GT.T) RETURN
      IFLAG = 2
      IF (E.GT.EPS.AND.E.LT.T) RETURN
      IFLAG = 1
      IF (E.LT.EPS.AND.E.GT.T) RETURN
      IFLAG = 0
      RETURN
      END SUBROUTINE Q1DAX

      SUBROUTINE DAMAX(N,SX,INCX,IDAMAX)
      DOUBLE PRECISION SX(*),SMAX,XMAG
      INTEGER N,INCX,I,II,NS,IDAMAX
      IDAMAX = 0
      IF(N.LE.0) RETURN
      IDAMAX = 1
      IF(N.LE.1)RETURN
      IF(INCX.EQ.1)GOTO 20
!     CODE FOR INCREMENTS NOT EQUAL TO 1.
      SMAX = DABS(SX(1))
      NS = N*INCX
      II = 1
          DO 10 I=1,NS,INCX
          XMAG = DABS(SX(I))
          IF(XMAG.LE.SMAX) GO TO 5
          IDAMAX = II
          SMAX = XMAG
    5     II = II + 1
   10     CONTINUE
      RETURN
!     CODE FOR INCREMENTS EQUAL TO 1.
   20 SMAX = DABS(SX(1))
      DO 30 I = 2,N 
         XMAG = DABS(SX(I))
         IF(XMAG.LE.SMAX) GO TO 30
         IDAMAX = I 
         SMAX = XMAG
   30 CONTINUE
      RETURN
      END SUBROUTINE DAMAX

      SUBROUTINE GL15T (F, A, B, XL, XR, R, AE, RA, RASC, FMIN, FMAX,   &
      RPAR,IPAR)

      IMPLICIT NONE
      DOUBLE PRECISION A, AE, B, DHLGTH, EPMACH, FC, FMAX, FMIN,        &
      FVAL1, FVAL2, HLGTH, PHI, PHIP, PHIU, R, RA, RASC, RESG,          &
      RESK, RESKH, SL, SR, UFLOW, WG, WGK, XGK, FSUM, FV1, FV2,         &
      XL, XR, CENTR, ABSC, U, RPAR
      INTEGER J, JTW, JTWM1, IPAR
      DIMENSION RPAR(1),IPAR(1)
      SAVE EPMACH, UFLOW
      DIMENSION FV1(7), FV2(7), WG(4), WGK(8), XGK(8)
      EXTERNAL F
!
!     THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1)
!     BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!     CORRESPONDING WEIGHTS ARE GIVEN.
!     XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
!              XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
!              GAUSS RULE
!              XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!              ADDED TO THE 7-POINT GAUSS RULE
!     WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
!     WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
!
      DATA WG (1) / 0.129484966168869693270611432679082D0 /
      DATA WG (2) / 0.279705391489276667901467771423780D0 /
      DATA WG (3) / 0.381830050505118944950369775488975D0 /
      DATA WG (4) / 0.417959183673469387755102040816327D0 /
!
      DATA XGK (1) / 0.991455371120812639206854697526329D0 /
      DATA XGK (2) / 0.949107912342758524526189684047851D0 /
      DATA XGK (3) / 0.864864423359769072789712788640926D0 /
      DATA XGK (4) / 0.741531185599394439863864773280788D0 /
      DATA XGK (5) / 0.586087235467691130294144838258730D0 /
      DATA XGK (6) / 0.405845151377397166906606412076961D0 /
      DATA XGK (7) / 0.207784955007898467600689403773245D0 /
      DATA XGK (8) / 0.000000000000000000000000000000000D0 /
!
      DATA WGK (1) / 0.022935322010529224963732008058970D0 /
      DATA WGK (2) / 0.063092092629978553290700663189204D0 /
      DATA WGK (3) / 0.104790010322250183839876322541518D0 /
      DATA WGK (4) / 0.140653259715525918745189590510238D0 /
      DATA WGK (5) / 0.169004726639267902826583426598550D0 /
      DATA WGK (6) / 0.190350578064785409913256402421014D0 /
      DATA WGK (7) / 0.204432940075298892414161999234649D0 /
      DATA WGK (8) / 0.209482141084727828012999174891714D0 /
!
      PHI(U) = XR - (XR - XL) * U * U * (2.D0 * U + 3.D0)
      PHIP(U) = - 6.D0 * U * (U + 1.D0)
!
!           LIST OF MAJOR VARIABLES:
!
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - R OF THE 7-POINT GAUSS FORMULA
!           RESK   - R OF THE 15-POINT KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!                    I.E. TO I/(B-A)
!
!           MACHINE DEPENDENT CONSTANTS:
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
      DATA EPMACH, UFLOW / 0.0D0, 0.0D0 /
      IF (EPMACH.EQ.0.0D0) THEN
         EPMACH = 2.22D-16
         UFLOW = 2.23D-308
      ENDIF
!
      IF (XL.LT.XR) THEN
         SL = SNGL(XL)
         SR = SNGL(XR)
      ELSE
         SL = SNGL(XR)
         SR = SNGL(XL)
      ENDIF
      HLGTH = 0.5D+00 * (B - A)
      CENTR = A + HLGTH
      DHLGTH = DABS(HLGTH)
!
!     COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
!     THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
      U = (CENTR - XR) / (XR - XL)
      PHIU = PHI(U)
      IF (PHIU.LE.SL.OR.PHIU.GE.SR) PHIU = CENTR
      CALL F(PHIU,RPAR,IPAR,FMIN)
      FMAX = FMIN
      FC = FMIN * PHIP(U)
      RESG = FC * WG(4)
      RESK = FC * WGK(8)
      RA = DABS(RESK)
      DO 10 J = 1, 3
         JTW = J * 2
         ABSC = HLGTH * XGK(JTW)
         U = (CENTR - ABSC - XR) / (XR - XL)
         PHIU = PHI(U)
         IF (PHIU.LE.SL.OR.PHIU.GE.SR) PHIU = CENTR
         CALL F(PHIU,RPAR,IPAR,FVAL1)
         FMAX = MAX(FMAX,FVAL1)
         FMIN = MIN(FMIN,FVAL1)
         FVAL1 = FVAL1 * PHIP(U)
         U = (CENTR + ABSC - XR) / (XR - XL)
         PHIU = PHI(U)
         IF (PHIU.LE.SL.OR.PHIU.GE.SR) PHIU = CENTR
         CALL F(PHIU,RPAR,IPAR,FVAL2)
         FMAX = MAX(FMAX,FVAL2)
         FMIN = MIN(FMIN,FVAL2)
         FVAL2 = FVAL2 * PHIP(U)
         FV1(JTW) = FVAL1
         FV2(JTW) = FVAL2
         FSUM = FVAL1 + FVAL2
         RESG = RESG + WG(J) * FSUM
         RESK = RESK + WGK(JTW) * FSUM
         RA = RA + WGK(JTW) * (DABS(FVAL1) + DABS(FVAL2) )
   10 END DO
      DO 15 J = 1, 4
         JTWM1 = J * 2 - 1
         ABSC = HLGTH * XGK(JTWM1)
         U = (CENTR - ABSC - XR) / (XR - XL)
         PHIU = PHI(U)
         IF (PHIU.LE.SL.OR.PHIU.GE.SR) PHIU = CENTR
         CALL F(PHIU,RPAR,IPAR,FVAL1)
         FMAX = MAX(FMAX,FVAL1)
         FMIN = MIN(FMIN,FVAL1)
         FVAL1 = FVAL1 * PHIP(U)
         U = (CENTR + ABSC - XR) / (XR - XL)
         PHIU = PHI(U)
         IF (PHIU.LE.SL.OR.PHIU.GE.SR) PHIU = CENTR
         CALL F(PHIU,RPAR,IPAR,FVAL2)
         FMAX = MAX(FMAX,FVAL2)
         FMIN = MIN(FMIN,FVAL2)
         FVAL2 = FVAL2 * PHIP(U)
         FV1(JTWM1) = FVAL1
         FV2(JTWM1) = FVAL2
         FSUM = FVAL1 + FVAL2
         RESK = RESK + WGK(JTWM1) * FSUM
         RA = RA + WGK(JTWM1) * (DABS(FVAL1) + DABS(FVAL2) )
   15 END DO
      RESKH = RESK * 0.5D+00
      RASC = WGK(8) * DABS(FC-RESKH)
      DO 20 J = 1, 7
         RASC = RASC + WGK(J) * (DABS(FV1(J) - RESKH) + DABS(FV2(J)&
         - RESKH))
   20 END DO
      R = RESK * HLGTH
      RA = RA * DHLGTH
      RASC = RASC * DHLGTH
      AE = DABS((RESK - RESG) * HLGTH)
      IF (RASC.NE.0.0D+00.AND.AE.NE.0.0D+00) AE = RASC * MIN(0.1D+01,  &
      (0.2D+03 * AE / RASC) **1.5D+00)
      IF (RA.GT.UFLOW / (0.5D+02 * EPMACH) ) AE = MAX ((EPMACH *       &
      0.5D+02) * RA, AE)
      RETURN
      END SUBROUTINE GL15T

END MODULE Q1DAMODULE
