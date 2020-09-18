      SUBROUTINE POWSAV(HESS, GRAD, XPARAM, PMAT, ILOOP, BMAT, IPOW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION HESS(97,*),GRAD(*),BMAT(97,*),IPOW(9), XPARAM(*)
     +, PMAT(*)
**********************************************************************
*
* POWSAV STORES AND RESTORES DATA USED IN THE SIGMA GEOMETRY
*        OPTIMISATION.
*
*  ON INPUT HESS   = HESSIAN MATRIX, PARTIAL OR WHOLE.
*           GRAD   = GRADIENTS.
*           XPARAM = CURRENT STATE OF PARAMETERS.
*           ILOOP  = INDEX OF HESSIAN, OR FLAG OF POINT REACHED SO-FAR.
*           BMAT   = "B" MATRIX!
*           IPOW   = INDICES AND FLAGS.
*           IPOW(9)= 0 FOR RESTORE, 1 FOR DUMP
*
**********************************************************************
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR), DUMY(MAXPAR)
      COMMON /DENSTY/ P(MPACK), PA(MPACK), PB(MPACK)
      COMMON /ALPARM/ ALPARM(3,MAXPAR),JLOOP,X0, X1, X2
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     +                NLAST(NUMATM), NORBS, NELECS,
     1                NALPHA, NBETA, NCLOSE, NOPEN
      COMMON /PATH  / LATOM,LPARAM,REACT(100)
      IR=9
      REWIND IR
      REWIND 10
      IF(IPOW(9) .EQ. 1) THEN
        WRITE(6,'(//10X,''- - - - - - - TIME UP - - - - - - -'',//)')
        WRITE(6,'(//10X,'' - THE CALCULATION IS BEING DUMPED TO DISK'',
     +  /10X,''   RESTART IT USING THE MAGIC WORD "RESTART"'')')
        FUNCT1=SQRT(DOT(GRAD,GRAD,NVAR))
        WRITE(6,'(//10X,''CURRENT VALUE OF GRADIENT NORM =''
     +  ,F12.6)')FUNCT1
        WRITE(6,'(/10X,''CURRENT VALUE OF GEOMETRY'',/)')
        CALL GEOUT
        WRITE(IR)IPOW,ILOOP
        WRITE(IR)(XPARAM(I),I=1,NVAR)
        WRITE(IR)(  GRAD(I),I=1,NVAR)
        WRITE(IR)((HESS(J,I),J=1,NVAR),I=1,NVAR)
        WRITE(IR)((BMAT(J,I),J=1,NVAR),I=1,NVAR)
        LINEAR=(NVAR*(NVAR+1))/2
        WRITE(IR)(PMAT(I),I=1,LINEAR)
        LINEAR=(NORBS*(NORBS+1))/2
        WRITE(10)(PA(I),I=1,LINEAR)
        IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
        IF(LATOM .NE. 0) THEN
          WRITE(IR)((ALPARM(J,I),J=1,3),I=1,NVAR)
          WRITE(IR)JLOOP,X0, X1, X2
        ENDIF
        STOP
      ELSE
        WRITE(6,'(//10X,'' RESTORING DATA FROM DISK''/)')
        READ(IR)IPOW,ILOOP
        READ(IR)(XPARAM(I),I=1,NVAR)
        READ(IR)(  GRAD(I),I=1,NVAR)
        READ(IR)((HESS(J,I),J=1,NVAR),I=1,NVAR)
        READ(IR)((BMAT(J,I),J=1,NVAR),I=1,NVAR)
        FUNCT1=SQRT(DOT(GRAD,GRAD,NVAR))
        WRITE(6,'(10X,''FUNCTION ='',F13.6//)')FUNCT1
        LINEAR=(NVAR*(NVAR+1))/2
        READ(IR)(PMAT(I),I=1,LINEAR)
        LINEAR=(NORBS*(NORBS+1))/2
        READ(10)(PA(I),I=1,LINEAR)
        IF(NALPHA.NE.0)READ(10)(PB(I),I=1,LINEAR)
        IF(LATOM.NE.0) THEN
          READ(IR)((ALPARM(J,I),J=1,3),I=1,NVAR)
          READ(IR)JLOOP,X0, X1, X2
          ILOOP=ILOOP+1
        ENDIF
        ILOOP=ILOOP+1
        RETURN
      ENDIF
      END