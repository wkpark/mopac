      SUBROUTINE POWSQ(XPARAM, NVAR, FUNCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XPARAM(*)
      COMMON /MESAGE/ IFLEPO,ISCF
**********************************************************************
*
*   POWSQ OPTIMISES THE GEOMETRY BY MINIMISING THE GRADIENT NORM.
*         THUS BOTH GROUND AND TRANSITION STATE GEOMETRIES CAN BE 
*         CALCULATED. IT IS ROUGHLY EQUIVALENT TO FLEPO, FLEPO MINIMISES
*         THE ENERGY, POWSQ MINIMISES THE GRADIENT NORM. 
*
*  ON ENTRY XPARAM = VALUES OF PARAMETERS TO BE OPTIMISED.
*           NVAR   = NUMBER OF PARAMETERS TO BE OPTIMISED.
*
*  ON EXIT  XPARAM = OPTIMISED PARAMETERS.
*           FUNCT  = HEAT OF FORMATION IN KCALS.
*
**********************************************************************
C        *****  ROUTINE PERFORMS  A LEAST SQUARES MINIMIZATION  *****
C        *****  OF A FUNCTION WHICH IS A SUM OF SQUARES.        *****
C        *****  INITIALLY WRITTEN BY J.W. MCIVER JR. AT SUNY/   *****
C        *****  BUFFALO, SUMMER 1971.  REWRITTEN AND MODIFIED   *****
C        *****  BY A.K. AT SUNY BUFFALO AND THE UNIVERSITY OF   *****
C        *****  TEXAS.  DECEMBER 1973                           *****
C
      COMMON /GEOVAR/ NDUM,LOC(2,MAXPAR),XARAM(MAXPAR)
      COMMON /GMETRY/ GEO(3,NUMATM)
      COMMON /KEYWRD/ KEYWRD
      COMMON /TIME  / TIME0
      COMMON /GEOSYM/ NDEP, LOCPAR(200), IDEPFN(200), LOCDEP(200)
      COMMON /GRADNT/ GRAD(MAXPAR),GNFINA
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     +                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     +                NLAST(NUMATM), NORBS, NELECS,
     1                NALPHA, NBETA, NCLOSE, NOPEN
      COMMON /NUMCAL/ NUMCAL
      COMMON /SIGMA1/ GNEXT, AMIN, ANEXT
      COMMON /SIGMA2/ GNEXT1(MAXPAR), GMIN1(MAXPAR)
      COMMON /NLLCOM/ HESS,BMAT,PMAT
      DIMENSION HESS(MAXPAR,MAXPAR), BMAT(MAXPAR,MAXPAR)
      DIMENSION GNORM(MAXPAR), IPOW(9), SIG(MAXPAR),
     1          PMAT(MAXHES), E1(MAXPAR), E2(MAXPAR),
     2          GREF1(MAXPAR), P(MAXPAR), WORK(MAXPAR), 
     3          PVEC(MAXPAR), EIG(MAXPAR), Q(MAXPAR)
      DIMENSION ISWAP(3), XOLD(MAXPAR)
      LOGICAL DEBUG, RESTRT, TIMES, OKC, OKF, ROUGH, REBEGN, NOPRT
      CHARACTER*80 KEYWRD
      DATA  ICALCN /0/, ISWAP /2,3,1/
      IF(ICALCN.NE.NUMCAL) THEN
          ICALCN=NUMCAL
          RESTRT=(INDEX(KEYWRD,'RESTART') .NE. 0)
          REBEGN=(INDEX(KEYWRD,'REBEGIN') .NE. 0)
          ROUGH=.FALSE.
          TIME1=SECOND()
          TIME2=TIME1
          TIMES=(INDEX(KEYWRD,'TIME') .NE. 0)
          TOTIME=3600
          I=INDEX(KEYWRD,' T=') 
          IF(I.NE.0) TOTIME=READA(KEYWRD,I)
          IF(I.NE.0)
     +    WRITE(6,'(//10X,'' TIME FOR THIS STEP ='',F8.2)')TOTIME
          STEP=0.02D0
          ILOOP=1
          NAT3=NUMAT*3
          XINC=0.00529167D0
          RHO2=1.D-8
          TOL2=2.D-1
          IF(INDEX(KEYWRD,'PREC') .NE. 0) TOL2=1.D-2
          IF(INDEX(KEYWRD,'VPREC') .NE. 0) TOL2=1.D-3
          DEBUG = (INDEX(KEYWRD,'POWSQ') .NE. 0)
          IF(RESTRT) THEN
C
C   RESTORE STORED DATA
C
              IPOW(9)=0
              CALL POWSAV(HESS,GMIN1,XPARAM,PMAT,ILOOP,BMAT,IPOW)
              DO 56 I=1,NVAR
              GRAD(I)=GMIN1(I)
  56          GNEXT1(I)=GMIN1(I)
              WRITE(6,'('' XPARAM'',6F12.6)')(XPARAM(I),I=1,NVAR)
              IF(REBEGN)ILOOP=1
              IF(ILOOP .GT. 0) THEN
               WRITE(6,'(//10X,'' RESTARTING AT POINT'',I3)')ILOOP
              ELSE
               WRITE(6,'(//10X,''RESTARTING IN OPTIMISATION'',
     +         '' ROUTINES'')')
              ENDIF
          ENDIF
*
*   DEFINITIONS:   NVAR   = NUMBER OF GEOMETRIC VARIABLES = 3*NUMAT-6
*
      ENDIF
      ILOOP=MIN(1,NVAR)
      NOPRT=(NVAR.LT.0)
      NVAR=ABS(NVAR)
      IF(DEBUG) THEN
          WRITE(6,'('' XPARAM'')')
      WRITE(6,'(5(2I3,F10.4))')(LOC(1,I),LOC(2,I),XPARAM(I),I=1,NVAR)
      ENDIF
      IF( .NOT. RESTRT .OR. REBEGN)
     +    CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
      IF(DEBUG) THEN
              WRITE(6,'('' STARTING GRADIENTS'')')
              WRITE(6,'(3X,8F9.4)')(GRAD(I),I=1,NVAR)
      ENDIF
      GMIN=SQRT(DOT(GRAD,GRAD,NVAR))
      GLAST=GMIN  
      DO 7 I=1,NVAR
      GNEXT1(I)=GRAD(I)
      GMIN1(I)=GNEXT1(I)
   7  CONTINUE
C
C    NOW TO CALCULATE THE HESSIAN MATRIX.
C
      IF(ILOOP.LT.0) GOTO 19
C
C   CHECK THAT HESSIAN HAS NOT ALREADY BEEN CALCULATED.
C
          DO 12 ILOOP=ILOOP,NVAR
          TIME1=SECOND()
          XPARAM(ILOOP)=XPARAM(ILOOP) + XINC
          CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
          IF(DEBUG)WRITE(6,'(I3,12(8F9.4,/3X))')
     +    ILOOP,(GRAD(IF),IF=1,NVAR)
          GRAD(ILOOP)=GRAD(ILOOP)+1.D-5
          XPARAM(ILOOP)=XPARAM(ILOOP) - XINC
          DO 1234 J=1,NVAR
 1234         HESS(ILOOP,J)=-(GRAD(J)-GNEXT1(J))/XINC
          TIME2=SECOND()
          TSTEP=TIME2-TIME1
          IF(TIMES)WRITE(6,'('' TIME FOR STEP:'',F8.2,'' LEFT'',F8.2)')
     +    TSTEP, TOTIME-TIME2+TIME0
          IF( TOTIME-TIME2+TIME0 .LT. TSTEP*2.D0) THEN
C
C  STORE RESULTS TO DATE.
C
              IPOW(9)=1
              CALL POWSAV(HESS,GMIN1,XPARAM,PMAT,ILOOP,BMAT,IPOW)
          ENDIF
  12      CONTINUE
      ILOOP=PARAM+1
C        *****  SCALE -HESSIAN- MATRIX                           *****
      IF( DEBUG) THEN
      WRITE(6,'(//10X,''UN-NORMALISED HESSIAN MATRIX'')')
      DO 132 I=1,NVAR
  132 WRITE(6,'(8F10.4)')(HESS(J,I),J=1,NVAR)
      ENDIF
      DO 14 I=1,NVAR
      SUM = 0.0D0
      DO 13 J=1,NVAR
  13  SUM = SUM+HESS(I,J)**2
  14  WORK(I) = 1.0D0/SQRT(SUM)
      DO 15 I=1,NVAR
      DO 15 J=1,NVAR
  15  HESS(I,J) = HESS(I,J)*WORK(I)
      IF( DEBUG) THEN
      WRITE(6,'(//10X,''HESSIAN MATRIX'')')
      DO 131 I=1,NVAR
  131 WRITE(6,'(8F10.4)')(HESS(J,I),J=1,NVAR)
      ENDIF
C        *****  INITIALIZE B MATIRX                        *****
      DO 16 I=1,NVAR
      DO 116 J=1,NVAR
  116 BMAT(I,J) = 0.0D0
  16  BMAT(I,I) = WORK(I)*2.D0
***************************************************************************
*
*  THIS IS THE START OF THE BIG LOOP TO OPTIMISE THE GEOMETRY
*
***************************************************************************
      ILOOP=-99
      TSTEP=TSTEP*4
  18  CONTINUE
      IF( TOTIME-TIME2+TIME0 .LT. TSTEP*2.D0) THEN
C
C  STORE RESULTS TO DATE.
C
          IPOW(9)=1
          CALL POWSAV(HESS,GMIN1,XPARAM,PMAT,ILOOP,BMAT,IPOW)
      ENDIF
  19  CONTINUE
C        *****  FORM-A- DAGGER-A- IN PA SLONG WITH -P-     *****
      IJ=0
      DO 20 J=1,NVAR
      DO 20 I=1,J
      IJ=IJ+1
      SUM = 0.0D0
      DO 17 K=1,NVAR
  17  SUM = SUM + HESS(I,K)*HESS(J,K)
  20  PMAT(IJ) = SUM
      DO 21 I=1,NVAR
      SUM = 0.0D0
      DO 121 K=1,NVAR
  121 SUM = SUM-HESS(I,K)*GMIN1(K)
  21  P(I) = -SUM
      L=0
      IF(DEBUG) THEN
      WRITE(6,'(/10X,''P MATRIX IN POWSQ'')')
      CALL VECPRT(PMAT,NVAR)
      ENDIF
      CALL HQRII(PMAT,NVAR,NVAR,EIG,PVEC)

C        *****  CHECK FOR ZERO EIGENVALUE                  *****
      IF(EIG(1).LT.RHO2) GO TO 25
      INDC = 2
C        *****  IF MATRIX IS NOT SINGULAR FORM INVERSE     *****
C        *****  BY BACK TRANSFORMING THE EIGENVECTORS      *****
      IJ=0
      DO 22 I=1,NVAR
      DO 22 J=1,I
      IJ=IJ+1
      SUM = 0.0D0
      DO 122 K=1,NVAR
  122 SUM = SUM+PVEC((K-1)*NVAR+J)*PVEC((K-1)*NVAR+I)/EIG(K)
  22  PMAT(IJ) = SUM
C        *****  FIND -Q- VECTOR                            *****
      L=0
      IL=L+1
      L=IL+I-1
      DO 24 I=1,NVAR
      SUM = 0.0D0
      DO 23 K=1,I
      IK=(I*(I-1))/2+K
  23  SUM = SUM+PMAT(IK)*P(K)
      IP1=I+1
      DO 26 K=IP1,NVAR
      IK=(K*(K-1))/2+I
  26  SUM=SUM+PMAT(IK)*P(K)
  24  Q(I) = SUM
      GO TO 31
  25  CONTINUE
C        *****  TAKE  -Q- VECTOR AS EIGENVECTOR OF ZERO     *****
C        *****  EIGENVALUE                                 *****
      DO 30 I=1,NVAR
  30  Q(I) = PVEC(I)
  31  CONTINUE
C        *****  FIND SEARCH DIRECTION                      *****
      DO 32 I=1,NVAR
      SIG(I) = 0.0D0
      DO 32 J=1,NVAR
  32  SIG(I) = SIG(I) + Q(J)*BMAT(I,J)
C        *****  DO A ONE DIMENSIONAL SEARCH                *****
      IF (DEBUG) THEN
      WRITE(6,'('' SEARCH VECTOR'')')
      WRITE(6,'(8F10.5)')(SIG(I),I=1,NVAR)
      ENDIF
C#      IF(ROUGH) THEN
C#          CALL LINMIN(XPARAM, ALPHA, SIG, NVAR, GMIN, OKC, OKF)
C#      ELSE
          CALL SEARCH(XPARAM, ALPHA, SIG, NVAR, GMIN, OKC, OKF)
          IF( NVAR .EQ. 1) GOTO 50
C
C  FIRST WE ATTEMPT TO OPTIMISE GEOMETRY USING SEARCH.
C  IF THIS DOES NOT WORK, THEN SWITCH TO LINMIN, WHICH ALWAYS WORKS,
C  BUT IS TWICE AS SLOW AS SEARCH.
C
          ROUGH=  (   .NOT.  OKF)
C#      ENDIF
      RMX = 0.0D0
      DO 33 K=1,NVAR
      RT = ABS(GMIN1(K))
      IF(RT.GT.RMX)RMX = RT
  33  CONTINUE
      IF(RMX.LT.TOL2) GO TO 50
C        *****  TWO STEP ESTIMATION OF DERIVATIVES         *****
      DO 34 K=1,NVAR
  34  E1(K) = (GMIN1(K)-GNEXT1(K))/(AMIN-ANEXT)
      RMU = DOT(E1,GMIN1,NVAR)/DOT(GMIN1,GMIN1,NVAR)
      DO 36 K=1,NVAR
  36  E2(K) = E1(K) - RMU*GMIN1(K)
C        *****  SCALE -E2- AND -SIG-                       *****
      SK = 1.0D0/SQRT(DOT(E2,E2,NVAR))
      DO 41 K=1,NVAR
  41  SIG(K) = SK*SIG(K)
      DO 42 K=1,NVAR
  42  E2(K) = SK*E2(K)
C        *****  FIND INDEX OF REPLACEMENT DIRECTION        *****
      PMAX = -1.0D+20
      DO 43 I=1,NVAR
      IF(ABS(P(I)*Q(I)).LE.PMAX) GO TO 43
      PMAX = ABS(P(I)*Q(I))
      ID = I
  43  CONTINUE
C        *****  REPLACE APPROPRIATE DIRECTION AND DERIVATIVE ***
      DO 44 K=1,NVAR
  44  HESS(ID,K) = -E2(K)
C        *****  REPLACE STARTING POINT                     *****
      DO 45 K=1,NVAR
  45  BMAT(K,ID) = SIG(K)/0.529167D0
      DO 46 K=1,NVAR
  46  GNEXT1(K) = GMIN1(K)
      IF( .NOT. NOPRT)WRITE(6,'('' GRADIENT ='',F15.6)')GMIN
      GLAST = GMIN
      IF( .NOT. NOPRT)WRITE(6,'(''  REPLACING DIRECTION '',I4)') ID
      INDC = 1
      TIME1=TIME2
      TIME2=SECOND()
      TSTEP=TIME2-TIME1
      IF(TIMES)WRITE(6,'('' TIME FOR STEP:'',F8.2,'' LEFT'',F8.2)')
     +TSTEP, TOTIME-TIME2+TIME0
      GO TO 18
  50  CONTINUE
      CALL COMPFG(XPARAM, .TRUE., FUNCT, .TRUE., GRAD, .TRUE.)
      DO 49 I=1,NVAR
  49  GRAD(I)=GMIN1(I)
      GNFINA=SQRT(DOT(GRAD,GRAD,NVAR))
      IFLEPO=11
      RETURN
      END
