      SUBROUTINE FORCE
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR),IDUMY,DUMY(MAXPAR)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),
     1                     LOCDEP(MAXPAR)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /FMATRX/ FMATRX(MAXPAR**2+MAXPAR*3+1),IDUMY2(4)
      COMMON /KEYWRD/ KEYWRD
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      PARAMETER (IPADD=2*MORB2+2*MAXORB-MAXPAR-MAXPAR*MAXPAR)
      COMMON /VECTOR/ CNORML(MAXPAR*MAXPAR),FREQ(MAXPAR),DUMMY(IPADD)
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /LAST  / LAST
      COMMON /MESAGE/ IFLEPO,ISCF
      COMMON /SYMOPS/ R(14,120), NSYM, IPO(NUMATM,120)
      COMMON /SIMBOL/ SIMBOL(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /COORD / COORD(3,NUMATM)
***********************************************************************
*
*   FORCE CALCULATES THE FORCE CONSTANTS FOR THE MOLECULE, AND THE
*         VIBRATIONAL FREQUENCIES.  ISOTOPIC SUBSTITUTION IS ALLOWED.
*
***********************************************************************
      COMMON /EULER / TVEC(3,3), ID
      COMMON /SCFTYP/ EMIN, LIMSCF
      COMMON /SCRACH/ STORE(MAXPAR**2)
      DIMENSION XPARAM(MAXPAR), GR(3,NUMATM),
     1DELDIP(3,MAXPAR), TRDIP(3,MAXPAR),LOCOLD(2,MAXPAR)
     2,REDMAS(MAXPAR), SHIFT(6), DIPT(MAXPAR), TRAVEL(MAXPAR)
     3, ROT(3,3), GEOREF(3,NUMATM), NAR(NUMATM), NBR(NUMATM),NCR(NUMATM)
      CHARACTER KEYWRD*241, KEYS(241)*1, ELEMNT*2, SIMBOL*10
      LOGICAL RESTRT, LINEAR, DEBUG, BARTEL, PRNT, LARGE, LIMSCF
      EQUIVALENCE (GRAD(1), GR(1,1)), (KEYWRD,KEYS(1))
C
C TEST GEOMETRY TO SEE IF IT IS OPTIMIZED
      TIME2=-1.D9
      CALL GMETRY(GEO,COORD)
      NVAOLD=NVAR
      DO 10 I=1,NVAR
         LOCOLD(1,I)=LOC(1,I)
   10 LOCOLD(2,I)=LOC(2,I)
      NVAR=0
      NDEOLD=NDEP
      NDEP=0
      NUMAT=0
      IF(LABELS(1) .NE. 99) NUMAT=1
      DO 30 I=2,NATOMS
         IF(LABELS(I).EQ.99) GOTO 30
         IF(I.EQ.2)ILIM=1
         IF(I.EQ.3)ILIM=2
         IF(I.GT.3)ILIM=3
C
C  IS IT A POLYMER?
C
         IF(LABELS(I).EQ.107) THEN
            ILIM=1
         ELSE
            NUMAT=NUMAT+1
         ENDIF
C$DOIT ASIS
         DO 20 J=1,ILIM
            NVAR=NVAR+1
            LOC(1,NVAR)=I
            LOC(2,NVAR)=J
   20    XPARAM(NVAR)=GEO(J,I)
   30 CONTINUE
C
C   IF A RESTART, THEN TSCF AND TDER WILL BE FAULTY, THEREFORE SET TO -1
C
      TSCF=-1.D0
      TDER=-1.D0
      PRNT=(INDEX(KEYWRD,'RC=') .EQ. 0)
      DEBUG=(INDEX(KEYWRD,'DFORCE') .NE. 0)
      LARGE=(INDEX(KEYWRD,'LARGE') .NE. 0)
      BARTEL=(INDEX(KEYWRD,'NLLSQ') .NE. 0)
      RESTRT=(INDEX(KEYWRD,'RESTART') .NE. 0)
      TIME1=SECOND()
      IF (RESTRT) THEN
C
C   CHECK TO SEE IF CALCULATION IS IN NLLSQ OR FORCE.
C
         IF(BARTEL)GOTO 50
C
C   CALCULATION IS IN FORCE
C
         GOTO 90
      ENDIF
      CALL COMPFG( XPARAM, .TRUE., ESCF, .TRUE., GRAD, .FALSE.)
      IF(PRNT)WRITE(6,'(//10X,''HEAT OF FORMATION ='',F12.6,
     1'' KCALS/MOLE'')')ESCF
      TIME2=SECOND()
      TSCF=TIME2-TIME1
      CALL COMPFG( XPARAM, .TRUE., ESCF1, .FALSE., GRAD, .TRUE.)
      TIME3=SECOND()
      TDER=TIME3-TIME2
      IF(PRNT)WRITE(6,'(//10X,''INTERNAL COORDINATE DERIVATIVES'',//3X,
     1''NUMBER  ATOM'',2X,''BOND'',9X,''  ANGLE'',10X,''DIHEDRAL'',/)')
      L=0
      IU=0
      DO 40 I=1,NATOMS
         IF(LABELS(I).EQ.99) GOTO 40
         L=L+1
         IL=IU+1
         IF(I .EQ. 1) IU=IL-1
         IF(I .EQ. 2) IU=IL
         IF(I .EQ. 3) IU=IL+1
         IF(I .GT. 3) IU=IL+2
         IF(LABELS(I).EQ.107)IU=IL
         IF(PRNT)WRITE(6,'(I6,4X,A2,F13.6,2F13.6)')
     1L,ELEMNT(LABELS(I)),(GRAD(J),J=IL,IU)
   40 CONTINUE
C   TEST SUM OF GRADIENTS
      GNORM=SQRT(DOT(GRAD,GRAD,NVAR))
      IF(PRNT)WRITE(6,'(//10X,''GRADIENT NORM ='',F10.5)') GNORM
      IF(GNORM.LT.10.D0) GOTO 70
      IF(INDEX(KEYWRD,' LET ') .NE. 0) THEN
         WRITE(6,'(///1X,''** GRADIENT IS VERY LARGE, BUT SINCE "LET"'',
     1'' IS USED, CALCULATION WILL CONTINUE'')')
         GOTO 90
      ENDIF
      WRITE(6,'(///1X,''** GRADIENT IS TOO LARGE TO ALLOW '',
     1    ''FORCE MATRIX TO BE CALCULATED, (LIMIT=10) **'',//)')
   50 CONTINUE
      DO 60 I=1,NVAR
   60 SIMBOL(I)='---'
      WRITE(6,'(//10X,'' GEOMETRY WILL BE OPTIMIZED FIRST'')')
      IF(BARTEL) THEN
         WRITE(6,'(15X,''USING NLLSQ'')')
         CALL NLLSQ(XPARAM,NVAR)
      ELSE
         WRITE(6,'(15X,''USING FLEPO'')')
         CALL FLEPO(XPARAM,NVAR,ESCF)
C
C  DID FLEPO USE ALL THE TIME ALLOWED?
C
         IF(IFLEPO.EQ.-1) RETURN
      ENDIF
      LIMSCF=.FALSE.
      CALL COMPFG( XPARAM, .TRUE., ESCF, .TRUE., GRAD, .TRUE.)
      CALL WRITMO(TIME1,ESCF)
      WRITE(6,'(//10X,''GRADIENT NORM ='',F10.7)') GNORM
      CALL GMETRY(GEO,COORD)
   70 CONTINUE
         DO 80 J=1,NATOMS
      NAR(J)=NA(J)
      NBR(J)=NB(J)
      NCR(J)=NC(J)
      DO 80 I=1,3
   80 GEOREF(I,J)=GEO(I,J)
C
C NOW TO CALCULATE THE FORCE MATRIX
C
C CHECK OUT SYMMETRY
   90 CONTINUE
C
C   NEED TO ENSURE THAT XYZINT WILL WORK CORRECTLY BEFORE CALL
C   TO DRC.
C
      L=0
      DO 100 I=1,NATOMS
         IF(LABELS(I).NE.99)THEN
            L=L+1
            LABELS(L)=LABELS(I)
         ENDIF
  100 CONTINUE
      NATOMS=NUMAT
      CALL XYZINT(COORD,NUMAT,NA,NB,NC,1.D0,GEO)
      CALL GMETRY(GEO,COORD)
      IF(INDEX(KEYWRD,'THERMO').NE.0 .AND.GNORM.GT.1.D0) THEN
         WRITE(6,'(//30X,''**** WARNING ****'',//
     110X,'' GRADIENT IS VERY LARGE FOR A THERMO CALCULATION'',/
     210X,'' RESULTS ARE LIKELY TO BE INACCURATE IF THERE ARE'')')
         WRITE(6,'(10X,'' ANY LOW-LYING VIBRATIONS (LESS THAN ABOUT ''
     1,''400CM-1)'')')
         WRITE(6,'(10X,'' GRADIENT NORM SHOULD BE LESS THAN ABOUT '',
     1''0.2 FOR THERMO'',/10X,'' TO GIVE ACCURATE RESULTS'')')
      ENDIF
      IF(TSCF.GT.0.D0) THEN
         WRITE(6,'(//10X,''TIME FOR SCF CALCULATION ='',F8.2)')TSCF
         WRITE(6,'(//10X,''TIME FOR DERIVATIVES     ='',F8.2)')TDER
      ENDIF
      IF(NDEP.GT.0) THEN
         WRITE(6,'(//10X,''SYMMETRY WAS SPECIFIED, BUT '',
     1''CANNOT BE USED HERE'')')
         NDEP=0
      ENDIF
      IF(PRNT)CALL AXIS(COORD,NUMAT,A,B,C,WTMOL,2,ROT)
      NVIB=3*NUMAT-6
      IF(ABS(C).LT.1.D-20)NVIB=NVIB+1
      IF(ID.NE.0)NVIB=3*NUMAT-3
      IF(PRNT) THEN
         WRITE(6,'(/9X,''ORIENTATION OF MOLECULE IN FORCE CALCULATION'')
     1')
         WRITE(6,'(/,4X,''NO.'',7X,''ATOM'',9X,''X'',
     19X,''Y'',9X,''Z'',/)')
      ENDIF
      L=0
      DO 110 I=1,NATOMS
         IF(LABELS(I) .EQ. 99) GOTO 110
         L=L+1
         IF(PRNT)WRITE(6,'(I6,7X,I3,4X,3F10.4)')
     1    L,LABELS(I),(COORD(J,L),J=1,3)
  110 CONTINUE
      CALL FMAT(FMATRX, NVIB, TSCF, TDER, DELDIP,ESCF)
      NA(1)=0
         DO 120 J=1,NATOMS
      NA(J)=NAR(J)
      NB(J)=NBR(J)
      NC(J)=NCR(J)
      DO 120 I=1,3
  120 GEO(I,J)=GEOREF(I,J)
      IF(NVIB.LT.0)THEN
         NDEP=NDEOLD
         NVAR=0
         RETURN
      ENDIF
C
C   THE FORCE MATRIX IS PRINTED AS AN ATOM-ATOM MATRIX RATHER THAN
C   AS A 3N*3N MATRIX, AS THE 3N MATRIX IS VERY CONFUSING!
C
      IJ=0
      IU=0
      DO 150 I=1,NUMAT
         IL=IU+1
         IU=IL+2
         IM1=I-1
         JU=0
         DO 140 J=1,IM1
            JL=JU+1
            JU=JL+2
            SUM=0.D0
C$DOIT ASIS
            DO 130 II=IL,IU
C$DOIT ASIS
               DO 130 JJ=JL,JU
  130       SUM=SUM+FMATRX((II*(II-1))/2+JJ)**2
            IJ=IJ+1
  140    STORE(IJ)=SQRT(SUM)
         IJ=IJ+1
  150 STORE(IJ)=SQRT(
     1FMATRX(((IL+0)*(IL+1))/2)**2+
     2FMATRX(((IL+1)*(IL+2))/2)**2+
     3FMATRX(((IL+2)*(IL+3))/2)**2+2.D0*(
     4FMATRX(((IL+1)*(IL+2))/2-1)**2+
     5FMATRX(((IL+2)*(IL+3))/2-2)**2+
     6FMATRX(((IL+2)*(IL+3))/2-1)**2))
      IF(DEBUG) THEN
         WRITE(6,'(//10X,'' FULL FORCE MATRIX, INVOKED BY "DFORCE"'')')
         I=-NVAR
         CALL VECPRT(FMATRX,I)
      ENDIF
      IF(PRNT)THEN
         WRITE(6,'(//10X,'' FORCE MATRIX IN MILLIDYNES/ANGSTROM'')')
         CALL VECPRT(STORE,NUMAT)
      ENDIF
      L=(NVAR*(NVAR+1))/2
      DO 160 I=1,L
  160 STORE(I)=FMATRX(I)
      IF(PRNT) CALL AXIS(COORD,NUMAT,A,B,C,SUM,0,ROT)
      IF(PRNT)WRITE(6,'(//10X,''HEAT OF FORMATION ='',F12.6,
     1'' KCALS/MOLE'')')ESCF
      IF(LARGE)THEN
         CALL FRAME(STORE,NUMAT,0, SHIFT)
         CALL RSP(STORE,NVAR,NVAR,FREQ,CNORML)
         DO 170 I=NVIB+1,NVAR
            J=(FREQ(I)+50.D0)*0.01D0
  170    FREQ(I)=FREQ(I)-J*100
         IF(PRNT)THEN
            WRITE(6,'(//10X,''TRIVIAL VIBRATIONS, SHOULD BE ZERO'')')
            WRITE(6,'(/, F9.4,''=TX'',F9.4,''=TY'',F9.4,''=TZ'',
     1             F9.4,''=RX'',F9.4,''=RY'',F9.4,''=RZ'')')
     2(FREQ(I),I=NVIB+1,NVAR)
            WRITE(6,'(//10X,''FORCE CONSTANTS IN MILLIDYNES/ANGSTROM''
     1,'' (= 10**5 DYNES/CM)'',/)')
            WRITE(6,'(8F10.5)')(FREQ(I),I=1,NVIB)
C CONVERT TO WEIGHTED FMAT
            WRITE(6,'(//10X,'' ASSOCIATED EIGENVECTORS'')')
            I=-NVAR
            CALL MATOUT(CNORML,FREQ,NVIB,I,NVAR)
         ENDIF
      ENDIF
      CALL FREQCY(FMATRX,FREQ,CNORML,REDMAS,TRAVEL,.TRUE.,DELDIP)
C
C  CALCULATE ZERO POINT ENERGY
C
C
C  THESE CONSTANTS TAKEN FROM HANDBOOK OF CHEMISTRY AND PHYSICS 62ND ED.
C   N AVOGADRO'S NUMBER = 6.022045*10**23
C   H PLANCK'S CONSTANT = 6.626176*10**(-34)JHZ
C   C SPEED OF LIGHT    = 2.99792458*10**10 CM/SEC
C   CONST=0.5*N*H*C/(1000*4.184)
      CONST=1.4295718D-3
      SUM=0.D0
      DO 180 I=1,NVAR
  180 SUM=SUM+FREQ(I)
      SUM=SUM*CONST
      IF(PRNT)
     1WRITE(6,'(//10X,'' ZERO POINT ENERGY''
     2, F12.3,'' KILOCALORIES PER MOLE'')')SUM
      SUMM=0.D0
      DO 230 I=1,NVAR
         SUM1=1.D-20
C$DOIT VBEST
         DO 190 J=1,NVAR
  190    SUM1=SUM1+CNORML(J+(I-1)*NVAR)**2
         SUM1=1.D0/SQRT(SUM1)
C$DOIT ASIS
         DO 200 K=1,3
  200    GRAD(K)=0.D0
C$DOIT ASIS
         DO 220 K=1,3
            SUM=0.D0
C$DOIT VBEST
            DO 210 J=1,NVAR
  210       SUM=SUM+CNORML(J+(I-1)*NVAR)*DELDIP(K,J)
            SUMM=SUMM+ABS(SUM)
  220    TRDIP(K,I)=SUM*SUM1
         DIPT(I)=SQRT(TRDIP(1,I)**2+TRDIP(2,I)**2+TRDIP(3,I)**2)
  230 CONTINUE
      IF(PRNT)THEN
         WRITE(6,'(//3X,'' THE LAST'',I2,'' VIBRATIONS ARE THE'',
     1'' TRANSLATION AND ROTATION MODES'')')NVAR-NVIB
         WRITE(6,'(3X,'' THE FIRST THREE OF THESE BEING TRANSLATIONS'',
     1'' IN X, Y, AND Z, RESPECTIVELY'')')
      ENDIF
      IF(PRNT.AND.LARGE)THEN
         WRITE(6,'(//10X,'' FREQUENCIES, REDUCED MASSES AND '',
     1''VIBRATIONAL DIPOLES''/)')
         NTO6=NVAR/6
         NREM6=NVAR-NTO6*6
         IINC1=-5
         IF (NTO6.LT.1) GO TO 250
         DO 240 I=1,NTO6
            WRITE (6,'(/)')
            IINC1=IINC1+6
            IINC2=IINC1+5
            WRITE (6,'(3X,''I'',10I10)') (J,J=IINC1,IINC2)
            WRITE (6,'('' FREQ(I)'',6F10.4,/)') (FREQ(J),J=IINC1,IINC2)
            WRITE (6,'('' MASS(I)'',6F10.5,/)') (REDMAS(J),J=IINC1,IINC2
     1)
            WRITE (6,'('' DIPX(I)'',6F10.5)') (TRDIP(1,J),J=IINC1,IINC2)
            WRITE (6,'('' DIPY(I)'',6F10.5)') (TRDIP(2,J),J=IINC1,IINC2)
            WRITE (6,'('' DIPZ(I)'',6F10.5,/)') (TRDIP(3,J),J=IINC1,IINC
     12)
            WRITE (6,'('' DIPT(I)'',6F10.5)')
     1   (DIPT(J),J=IINC1,IINC2)
  240    CONTINUE
  250    CONTINUE
         IF (NREM6.LT.1) GO TO 260
         WRITE (6,'(/)')
         IINC1=IINC1+6
         IINC2=IINC1+(NREM6-1)
         WRITE (6,'(3X,''I'',10I10)') (J,J=IINC1,IINC2)
         WRITE (6,'('' FREQ(I)'',6F10.4)') (FREQ(J),J=IINC1,IINC2)
         WRITE (6,'(/,'' MASS(I)'',6F10.5)') (REDMAS(J),J=IINC1,IINC2)
         WRITE (6,'(/,'' DIPX(I)'',6F10.5)') (TRDIP(1,J),J=IINC1,IINC2)
         WRITE (6,'('' DIPY(I)'',6F10.5)') (TRDIP(2,J),J=IINC1,IINC2)
         WRITE (6,'('' DIPZ(I)'',6F10.5)') (TRDIP(3,J),J=IINC1,IINC2)
         WRITE (6,'(/,'' DIPT(I)'',6F10.5)')
     1   (DIPT(J),J=IINC1,IINC2)
  260    CONTINUE
      ENDIF
      IF(PRNT)THEN
         WRITE(6,'(//10X,'' NORMAL COORDINATE ANALYSIS'')')
         I=-NVAR
         CALL MATOUT(CNORML,FREQ,NVAR,I,NVAR)
      ENDIF
C
C   CARRY OUT IRC IF REQUESTED.
C
      IF(INDEX(KEYWRD,'IRC')+INDEX(KEYWRD,'DRC').eq.677)THEN
         DO 270 I=1,NVAR
            LOC(1,I)=0
  270    LOC(2,I)=0
         NVAR=NVAOLD
         DO 280 I=1,NVAR
            LOC(1,I)=LOCOLD(1,I)
  280    LOC(2,I)=LOCOLD(2,I)
         CALL XYZINT(COORD,NUMAT,NA,NB,NC,1.D0,GEO)
         LAST=1
         CALL DRC(CNORML,FREQ)
         NA(1)=0
         NDEP=NDEOLD
         NVAR=0
         DO 290 I=1,3
            DO 290 J=1,NATOMS
  290    GEO(I,J)=GEOREF(I,J)
         RETURN
      ENDIF
      CALL FREQCY(FMATRX,FREQ,CNORML,DELDIP,DELDIP,.FALSE.,DELDIP)
      WRITE(6,'(//10X,'' MASS-WEIGHTED COORDINATE ANALYSIS'')')
      I=-NVAR
      CALL MATOUT(CNORML,FREQ,NVAR,I,NVAR)
      CALL ANAVIB(COORD,FREQ,DIPT,NVAR,CNORML,STORE,
     1FMATRX,TRAVEL,REDMAS)
      IF(INDEX(KEYWRD,'THERMO').NE.0) THEN
         CALL GMETRY(GEO,COORD)
         I=INDEX(KEYWRD,' ROT')
         IF(I.NE.0) THEN
            SYM=READA(KEYWRD,I)
         ELSE
            SYM=1
         ENDIF
         LINEAR=(ABS(A*B*C) .LT. 1.D-10)
         I=INDEX(KEYWRD,' TRANS')
C
C   "I" IS GOING TO MARK THE BEGINNING OF THE GENUINE VIBRATIONS.
C
         IF(I.NE.0)THEN
            I=INDEX(KEYWRD,' TRANS=')
            IF(I.NE.0)THEN
               I=1+READA(KEYWRD,I)
               J=NVIB-I+1
               WRITE(6,'(//1X,''THE LOWEST'',I3,'' VIBRATIONS ARE NOT'',
     1/,'' TO BE USED IN THE THERMO CALCULATION'')')I-1
            ELSE
               WRITE(6,'(//10X,''SYSTEM IS A TRANSITION STATE'')')
               I=2
               J=NVIB-1
            ENDIF
         ELSE
            WRITE(6,'(//10X,''SYSTEM IS A GROUND STATE'')')
            I=1
            J=NVIB
         ENDIF
         CALL THERMO(A,B,C,LINEAR,SYM,WTMOL,FREQ(I),J,ESCF)
      ENDIF
      NA(1)=0
      NVAR=0
      NDEP=NDEOLD
      DO 300 I=1,3
         DO 300 J=1,NATOMS
  300 GEO(I,J)=GEOREF(I,J)
      RETURN
      END
