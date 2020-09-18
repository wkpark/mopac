      SUBROUTINE PATHS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON /PATH  / LATOM,LPARAM,REACT(100)
      COMMON /GEOVAR/ NVAR, LOC(2,MAXPAR), XPARAM(MAXPAR)
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /ALPARM/ ALPARM(3,MAXPAR),ILOOP, X0, X1, X2
*************************************************************************** 
*
*   PATH FOLLOWS A REACTION COORDINATE.   THE REACTION COORDINATE IS ON
*        ATOM LATOM, AND IS A DISTANCE IF LPARAM=1,
*                           AN ANGLE   IF LPARAM=2,
*                           AN DIHEDRALIF LPARAM=3.
*
*************************************************************************** 
      DIMENSION GD(MAXPAR),XLAST(MAXPAR),MDFP(20),XDFP(20)
      CHARACTER*80 KEYWRD
      CHARACTER*10 TYPE(3)
      DATA TYPE / 'ANGSTROMS ','DEGREES   ','DEGREES   '/
      ILOOP=1
      IF(INDEX(KEYWRD,'RESTAR') .NE. 0) THEN
          MDFP(9)=0
          CALL DFPSAV(TOTIME,XPARAM,GD,XLAST,FUNCT1,MDFP,XDFP)
          WRITE(6,'(//10X,'' RESTARTING AT POINT'',I3)')ILOOP
      ENDIF
          IF(ILOOP.GT.1) GOTO 20
          WRITE(6,'(''  ABOUT TO ENTER FLEPO FROM PATH'')')
          TIME1=SECOND()
          CALL FLEPO(XPARAM,NVAR,FUNCT)
      WRITE(6,'(''  OPTIMISED VALUES OF PARAMETERS, INITIAL POINT'')')
          CALL WRITE(TIME1,FUNCT)
  20      CONTINUE
          IF(ILOOP.GT.2) GOTO 30
          IF(ILOOP.EQ.1) THEN
              X0=REACT(1)
              X1=X0
              X2=REACT(2)
              IF(X2.LT. -100.D0) STOP
              GEO(LPARAM,LATOM)=X2
              DO 1 I=1,NVAR
              ALPARM(2,I)=XPARAM(I)
   1          ALPARM(1,I)=XPARAM(I)
              ILOOP=2
          ENDIF
          CALL FLEPO(XPARAM,NVAR,FUNCT)
      RNORD=REACT(2)
      IF(LPARAM.GT.1) RNORD=RNORD*57.29577951D0
      WRITE(6,'(1X,16(''*****'')//17X,''REACTION COORDINATE = ''
     +,F12.4,2X,A10,19X//1X,16(''*****''))')RNORD,TYPE(LPARAM)
          CALL WRITE(TIME1,FUNCT)
          TIME1=SECOND()
          DO 2 I=1,NVAR
   2      ALPARM(3,I)=XPARAM(I)
C
C   NOW FOR THE MAIN INTERPOLATION ROUTE
C
      IF(ILOOP.EQ.2)ILOOP=3
  30  CONTINUE
      DO 100 ILOOP = ILOOP,100
C
      IF(REACT(ILOOP).LT. -100.D0) THEN
           CALL EXIT
C
      ENDIF
C
      RNORD=REACT(ILOOP)
      IF(LPARAM.GT.1) RNORD=RNORD*57.29577951D0
      WRITE(6,'(1X,16(''*****'')//19X,''REACTION COORDINATE = ''
     +,F8.5,2X,A10,19X//1X,16(''*****''))')RNORD,TYPE(LPARAM)
C
      X3=REACT(ILOOP)
      C3=(X0**2-X1**2)*(X1-X2)-(X1**2-X2**2)*(X0-X1)
C      WRITE(6,'(''   C3:'',F13.7)')C3
      IF (ABS(C3) .LT. 1.D-8) THEN
C
C    WE USE A LINEAR INTERPOLATION
C
          CC1=0.D0
          CC2=0.D0
      ELSE
C    WE DO A QUADRATIC INTERPOLATION
C
          CC1=(X1-X2)/C3
          CC2=(X0-X1)/C3
      END IF
      CB1=1.D0/(X1-X2)
      CB2=(X1**2-X2**2)*CB1
C      WRITE(6,'(''   CB1'',F13.7)')CB1
C      WRITE(6,'(''   CB2'',F13.7)')CB2
C
C    NOW TO CALCULATE THE INTERPOLATED COORDINATES
C
      DO 4 I=1,NVAR
      DELF0=ALPARM(1,I)-ALPARM(2,I)
      DELF1=ALPARM(2,I)-ALPARM(3,I)
C      WRITE(6,'('' DELF0'',F13.7)') DELF0
C      WRITE(6,'('' DELF1'',F13.7)') DELF1
      ACONST = CC1*DELF0-CC2*DELF1
      BCONST = CB1*DELF1-ACONST*CB2
      CCONST = ALPARM(3,I) - BCONST*X2 - ACONST*X2**2
C      WRITE(6,'(''  ACONST,BCONST, CCONST'',3F13.7)')
C     + ACONST,BCONST,CCONST
      XPARAM(I)=CCONST+BCONST*X3+ACONST*X3**2
      ALPARM(1,I)=ALPARM(2,I)
   4  ALPARM(2,I)=ALPARM(3,I)
C      WRITE(6,'('' GUESSED PARAMETERS'')')
C      WRITE(6,'(8F10.6)')(XPARAM(I),I=1,NVAR)
      X0=X1
      X1=X2
      X2=X3
      GEO(LPARAM,LATOM)=REACT(ILOOP)
      CALL FLEPO(XPARAM,NVAR,FUNCT)
      CALL WRITE(TIME1,FUNCT)
      TIME1=SECOND()
      DO 5 I=1,NVAR
   5  ALPARM(3,I)=XPARAM(I)
  100 CONTINUE
      END
      
