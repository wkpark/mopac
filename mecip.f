      SUBROUTINE MECIP(COEFFS,NORBS,DELTAP, DELTA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION COEFFS(NORBS,NORBS), DELTAP(NMOS,NMOS),
     1 DELTA(NORBS,NMOS)
************************************************************************
*
*   MECIP WILL CORRECT THE TOTAL DENSITY MATRIX FOR THE EFFECT OF THE
*   C.I.
*              ON INPUT
*
*  COEFFS       : ALL M.O.'S (NORBS M.O.S)
*  NORBS        : NUMBER OF MOLECULAR ORBITALS = NUMBER OF A.O.'S
*  P            : TOTAL DENSITY MATRIX
*  NMOS         : NUMBER OF M.O.'S IN ACTIVE SPACE
*  VECTCI       : STATE VECTOR OF LENGTH LAB
*  MICROA(I,J)  : ALPHA OCCUPANCY OF M.O. 'I' IN MICROSTATE 'J'
*  MICROB(I,J)  : BETA  OCCUPANCY OF M.O. 'I' IN MICROSTATE 'J'
*
*  NOTE: THIS IS A MODIFICATION OF CODE ORIGINALLY WRITTEN BY
*        PROF. DANIEL LIOTARD
************************************************************************
      COMMON /CIBITS/ NMOS,LAB,NELEC, NBO(3)
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /NALMAT/ NALPHA(NMECI**2)
      COMMON /BASEOC/ OCCA(NMECI)
      COMMON /CIVECT/ VECTCI(NMECI**2),CONF(NMECI**4+1)
      COMMON /MICROS/ MICROA(NMECI,4*NMECI**2), MICROB(NMECI,4*NMECI**2)
C     INITIALIZE WITH THE OPPOSITE OF THE 'SCF' DENSITY.
      DO 10 I=1,NMOS
         DELTAP(I,I)=-OCCA(I)*2.D0
         DO 10 J=1,I-1
   10 DELTAP(I,J)=0.D0
C
C     ADD THE C.I. CORRECTION
      DO 120 ID=1,LAB
         DO 120 JD=1,ID
C     CHECK SPIN AGREEMENT
            IF(NALPHA(ID).NE.NALPHA(JD)) GO TO 120
            IX=0
            IY=0
            DO 20 J=1,NMOS
               IX=IX+ABS(MICROA(J,ID)-MICROA(J,JD))
   20       IY=IY+ABS(MICROB(J,ID)-MICROB(J,JD))
C     CHECK NUMBER OF DIFFERING M.O.
            IF(IX+IY.GT.2) GO TO 120
            IF(IX.EQ.2) THEN
C        DETERMINANTS ID AND JD DIFFER BY M.O I IN ID AND M.O J IN JD:
               DO 30 I=1,NMOS
   30          IF(MICROA(I,ID).NE.MICROA(I,JD)) GO TO 40
   40          IJ=MICROB(I,ID)
               DO 50 J=I+1,NMOS
                  IF(MICROA(J,ID).NE.MICROA(J,JD)) GO TO 60
   50          IJ=IJ+MICROA(J,ID)+MICROB(J,ID)
C        IJ GIVES THE SIGN OF THE PERMUTATION
   60          DELTAP(J,I)=DELTAP(J,I)+VECTCI(ID)*VECTCI(JD)*FLOAT(1-2*M
     1OD(IJ,2))
            ELSE IF(IY.EQ.2) THEN
C        DETERMINANTS ID AND JD DIFFER BY M.O J IN ID AND M.O I IN JD:
               DO 70 I=1,NMOS
   70          IF(MICROB(I,ID).NE.MICROB(I,JD)) GO TO 80
   80          IJ=0
               DO 90 J=I+1,NMOS
                  IF(MICROB(J,ID).NE.MICROB(J,JD)) GO TO 100
   90          IJ=IJ+MICROA(J,ID)+MICROB(J,ID)
  100          IJ=IJ+MICROA(J,ID)
               DELTAP(J,I)=DELTAP(J,I)+VECTCI(ID)*VECTCI(JD)*FLOAT(1-2*M
     1OD(IJ,2))
            ELSE
C        DETERMINANTS ID AND JD ARE IDENTICAL:
               DO 110 I=1,NMOS
  110          DELTAP(I,I)=DELTAP(I,I)+(MICROA(I,ID)+MICROB(I,ID))*VECTC
     1I(ID)**2
            ENDIF
  120 CONTINUE
C
C     BACK TRANSFORM INTO A.O. BASIS.
C     -------------------------------
C     P(C.I.) = P(SCF) + C * DELTAP * C'
      DO 130 I=1,NMOS
CDIR$ IVDEP
         DO 130 J=1,I-1
  130 DELTAP(J,I)=DELTAP(I,J)
C     STEP 1: DELTAP = C * DELTAP
      CALL MXM (COEFFS(1,NELEC+1),NORBS,DELTAP,NMOS,DELTA,NMOS)
C     STEP 2: P = P + DELTAP * C'
      IJ=0
      DO 150 I=1,NORBS
         DO 150 J=1,I
            IJ=IJ+1
            SUM=0.D0
            DO 140 K=1,NMOS
  140       SUM=SUM+DELTA(I,K)*COEFFS(J,NELEC+K)
  150 P(IJ)=P(IJ)+SUM
C     NOTE FROM D.L.: AT THIS POINT THE 'NATURAL ORBITALS' OF THIS STATE
C     CAN BE OBTAINED STRAIGHTWAY AS EIGENVECTORS OF THE DENSITY MATRIX.
      RETURN
      END
