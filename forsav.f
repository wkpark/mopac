      SUBROUTINE FORSAV(TIME,DELDIP,IPT,N3,FMATRX, COORD,NVAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES/NOLIST'
      DIMENSION FMATRX(*), DELDIP(3,*), COORD(*)
**************************************************************************
*
*  FORSAV SAVES AND RESTORES DATA USED IN THE FORCE CALCULATION.
*
* ON INPUT TIME   = TOTAL TIME ELAPSED SINCE THE START OF THE CALCULATION.
*          IPT    = LINE OF FORCE MATRIX REACHED, IF IN WRITE MODE,
*                 = 0 IF IN READ MODE.
*          FMATRX = FORCE MATRIX
**************************************************************************
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     +                NLAST(NUMATM), NORBS, NELECS,
     1                NALPHA, NBETA, NCLOSE, NOPEN
      IR=9
      IW=9
      REWIND 10
      REWIND IR
      IF( IPT .EQ. 0 ) THEN
C
C   READ IN FORCE DATA
C
      READ(IR)TIME,IPT
      LINEAR=(NVAR*(NVAR+1))/2
      READ(IR)(COORD(I),I=1,MAXPAR)
      READ(IR)(FMATRX(I),I=1,LINEAR)
      READ(IR)((DELDIP(J,I),J=1,3),I=1,IPT)
      LINEAR=(NORBS*(NORBS+1))/2
      READ(10)(PA(I),I=1,LINEAR)
      IF(NALPHA.NE.0)READ(10)(PB(I),I=1,LINEAR)
      RETURN
      ELSE
C
C    WRITE FORCE DATA
C
      REWIND IW
      WRITE(IW)TIME,IPT
      LINEAR=(NVAR*(NVAR+1))/2
      WRITE(IW)(COORD(I),I=1,MAXPAR)
      WRITE(IW)(FMATRX(I),I=1,LINEAR)
      WRITE(IW)((DELDIP(J,I),J=1,3),I=1,IPT)
      LINEAR=(NORBS*(NORBS+1))/2
      WRITE(10)(PA(I),I=1,LINEAR)
      IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
      IF(IPT.EQ.N3) THEN
          WRITE(6,'(//10X,''FORCE MATRIX WRITTEN TO DISK'')')
          ELSE
          STOP
      ENDIF
      ENDIF
      END
