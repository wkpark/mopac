      SUBROUTINE SEARCH(XPARAM,ALPHA,SIG,NVAR,GMIN,OKC,OKF, FUNCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XPARAM(*), SIG(*)
************************************************************************
*
* SEARCH PERFORMS A LINE SEARCH FOR POWSQ. IT MINIMIZES THE NORM OF
*        THE GRADIENT VECTOR IN THE DIRECTION SIG.
*
* ON INPUT  XPARAM = CURRENT POINT IN NVAR DIMENSIONAL SPACE.
*           ALPHA  = STEP SIZE (IN FACT ALPHA IS CALCULATED IN SEARCH).
*           SIG    = SEARCH DIRECTION VECTOR.
*           NVAR   = NUMBER OF PARAMETERS IN SIG (& XPARAM)
*
* ON OUTPUT XPARAM = PARAMETERS OF MINIMUM.
*           ALPHA  = DISTANCE TO MINIMUM.
*           GMIN   = GRADIENT NORM AT MINIMUM.
*           OKC    = EXITED BEFORE COUNTS WERE EXCEEDED (LOGICAL)
*           OKF    = FUNCTION WAS IMPROVED.
************************************************************************
      COMMON /SIGMA1/ GNEXT, AMIN, ANEXT
      COMMON /SIGMA2/  GNEXT1(MAXPAR), GMIN1(MAXPAR)
      COMMON/KEYWRD/ KEYWRD
      DIMENSION GRAD(MAXPAR),XREF(MAXPAR), GREF(MAXPAR), XMIN1(MAXPAR)
      CHARACTER*80 KEYWRD
      LOGICAL FIRST, DEBUG, OKC, OKF, NOPR
      DATA  FIRST /.TRUE./
      IF( FIRST ) THEN
         FIRST = .FALSE.
C
C    TOLG   = CRITERION FOR EXIT BY RELATIVE CHANGE IN GRADIENT.
C
         DEBUG=(INDEX(KEYWRD,'SEARCH') .NE. 0)
         NOPR=( .NOT. DEBUG)
         LOOKS=0
         OKF=.TRUE.
         TINY=0.1D0
         TOLERG=0.02D0
         IF(NVAR .EQ. 1)TIMY=0.D0
         G=100.D0
         XMAXM=2.D0
         ALPHA=0.1D0
      ENDIF
      ANEXT1=0.D0
      DO 10 I=1,NVAR
         GREF(I)  =GMIN1(I)
         GNEXT1(I)=GMIN1(I)
         XMIN1(I) =XPARAM(I)
   10 XREF(I)  =XPARAM(I)
      IF(ABS(ALPHA) .GT. 0.2)ALPHA=SIGN(0.2D0,ALPHA)
      IF(DEBUG) THEN
         WRITE(6,'('' SEARCH DIRECTION VECTOR'')')
         WRITE(6,'(6F12.6)')(SIG(I),I=1,NVAR)
         WRITE(6,'('' INITIAL GRADIENT VECTOR'')')
         WRITE(6,'(6F12.6)')(GMIN1(I),I=1,NVAR)
      ENDIF
      GB=DOT(GMIN1,GREF,NVAR)
      IF(DEBUG) WRITE(6,'('' GRADIENT AT START OF SEARCH:'',F16.6)')
     1SQRT(GB)
      GSTORE=GB
      AMIN=0.D0
      GMINN=1.D9
C
C
      TA=0.D0
      GA=GB
      GB=1.D9
      ITRYS=0
      GOTO 30
   20 SUM=GA/(GA-GB)
      ITRYS=ITRYS+1
      IF(ABS(SUM) .GT. 3.D0) SUM=SIGN(3.D0,SUM)
      ALPHA=(TB-TA)*SUM+TA
C
C         XPARAM IS THE GEOMETRY OF THE PREDICTED MINIMUM ALONG THE LINE
C
   30 CONTINUE
      DO 40 I=1,NVAR
   40 XPARAM(I)=XREF(I)+ALPHA*SIG(I)
C
C         CALCULATE GRADIENT NORM AND GRADIENTS AT THE PREDICTED MINIMUM
C
      DO 50 I=1,NVAR
   50 GRAD(I)=0.D0
      CALL COMPFG (XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
      LOOKS=LOOKS+1
C
C          G IS THE PROJECTION OF THE GRADIENT ALONG SIG.
C
      G=DOT(GREF,GRAD,NVAR)
      GTOT=SQRT(DOT(GRAD,GRAD,NVAR))
      IF( .NOT. NOPR)
     1WRITE(6,'('' LOOKS'',I3,'' ALPHA ='',F12.6,'' GRADIENT'',F12.3,
     2'' G  ='',F16.6)')
     3LOOKS,ALPHA,SQRT(DOT(GRAD,GRAD,NVAR)),G
      IF(GTOT .LT. GMINN) THEN
         GMINN=GTOT
         IF(ABS(AMIN-ALPHA) .GT.1.D-2) THEN
*
* WE CAN MOVE ANEXT TO A POINT NEAR, BUT NOT TOO NEAR, AMIN, SO THAT THE
* SECOND DERIVATIVESWILLBEREALISTIC(D2E/DX2=(GNEXT1-GMIN1)/(ANEXT-AMIN))
*
            ANEXT=AMIN
            DO 60 I=1,NVAR
   60       GNEXT1(I)=GMIN1(I)
         ENDIF
         AMIN=ALPHA
         DO 70 I=1,NVAR
            IF(GMINN.LT.GMIN) XMIN1(I)=XPARAM(I)
   70    GMIN1(I)=GRAD(I)
         IF(GMIN.GT.GMINN)GMIN=GMINN
      ENDIF
      IF(ITRYS .GT. 8) GOTO 80
      IF (ABS(G/GSTORE).LT.TINY .OR. ABS(G) .LT. TOLERG) GO TO 80
      IF(ABS(G) .LT. MAX(ABS(GA),ABS(GB)) .OR.
     1     GA*GB .GT. 0.D0 .AND. G*GA .LT. 0.D0) THEN
C
C   G IS AN IMPROVEMENT ON GA OR GB.
C
         IF(ABS(GB) .LT. ABS(GA))THEN
            TA=ALPHA
            GA=G
            GO TO 20
         ELSE
            TB=ALPHA
            GB=G
            GO TO 20
         ENDIF
      ELSE
C#         WRITE(6,'(//10X,'' FAILED IN SEARCH, SEARCH CONTINUING'')')
         GOTO 80
      ENDIF
   80 LNSTOP=4
      GMINN=SQRT(DOT(GMIN1,GMIN1,NVAR))
      DO 90 I=1,NVAR
   90 XPARAM(I)=XMIN1(I)
      IF(DEBUG) THEN
         WRITE(6,'('' AT EXIT FROM SEARCH'')')
         WRITE(6,'('' XPARAM'',6F12.6)')(XPARAM(I),I=1,NVAR)
         WRITE(6,'('' GNEXT1'',6F12.6)')(GNEXT1(I),I=1,NVAR)
         WRITE(6,'('' GMIN1 '',6F12.6)')(GMIN1(I),I=1,NVAR)
         WRITE(6,'('' AMIN, ANEXT,GMIN'',4F12.6)')
     1    AMIN,ANEXT,GMIN
      ENDIF
      RETURN
C
      END
