      SUBROUTINE WRITMO(TIME0,FUNCT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      CHARACTER KEYWRD*241
      DOUBLE PRECISION MECI
      COMMON /KEYWRD/ KEYWRD
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /HMATRX/ H(MPACK)
      COMMON /FOKMAT/ F(MPACK), FB(MPACK)
      COMMON /VECTOR/ C(MORB2),EIGS(MAXORB),CBETA(MORB2),EIGB(MAXORB)
      COMMON /DENSTY/ P(MPACK),PA(MPACK),PB(MPACK)
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                    LOCDEP(MAXPAR)
      COMMON / EULER/ TVEC(3,3), ID
      COMMON /RJKS  / RJKAB(NMECI,NMECI), RJKAA(NMECI,NMECI)
      COMMON /ERRFN / ERRFN(MAXPAR), AICORR(MAXPAR)
      COMMON /WORK1 /  FMAT2D(NPULAY*4), SEC(NPULAY*2), VEC(NPULAY*2),
     1                ALBAND(NPULAY*13)
      COMMON /PATH  / LATOM,LPARAM,REACT(200)
      COMMON /NUMCAL/ NUMCAL
      COMMON /NUMSCF/ NSCF
      COMMON /WMATRX/ WJ(N2ELEC), WK(N2ELEC)
      COMMON /ATHEAT/ ATHEAT
      COMMON /CORE  / CORE(107)
      COMMON /LAST  / LAST
      COMMON /SCRACH/ RXYZ(MPACK), XDUMY(MAXPAR**2-MPACK)
      COMMON /CIMATS/ ENGYCI(3),VECTCI(9),ECI(6)
      COMMON /MESAGE/ IFLEPO,IITER
      COMMON /ATMASS/ ATMASS(NUMATM)
      COMMON /ENUCLR/ ENUCLR
      COMMON /ELECT / ELECT
      COMMON /XYZGRA/ DXYZ(9*NUMATM)
      COMMON /GRADNT/ GRAD(MAXPAR), GNORM
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /GEOVAR/ NVAR, LOC(2,MAXPAR), IDUMY, XPARAM(MAXPAR)
      COMMON /BORN  / BP(NUMATM),FGB(NPACK),CCT1,ZEFF(NUMATM),
     1                QEFF(NUMATM)
      COMMON /SURF  / SURFCT,SURFACT(NUMATM),ATAR(NUMATM),ITYPE(NUMATM)
************************************************************************
*
*   WRITE PRINTS OUT MOST OF THE RESULTS.
*         IT SHOULD NOT ALTER ANY PARAMETERS, SO THAT IT CAN BE CALLED
*         AT ANY CONVENIENT TIME.
*
************************************************************************
      DIMENSION Q(MAXORB), Q2(MAXORB), COORD(3,NUMATM)
     1,IEL1(107), NELEMT(107), IEL2(107)
      DIMENSION W(N2ELEC), DUMY(3), EBP(NPACK), EST(100)
      LOGICAL UHF, CI, SINGLT, TRIPLT, EXCITD, PRTGRA, STILL, AQCHK, FAIL
      CHARACTER TYPE(3)*11, IDATE*24, CALCN(2)*5, GTYPE*13, GRTYPE*14,
     1          FLEPO(16)*58, ITER(2)*58, NUMBRS(11)*1, GETNAM*80
      CHARACTER*2 ELEMNT, IELEMT(20), CALTYP*7, NAMFIL*80
      SAVE ICALCN, NUMBRS, CALCN, TYPE, FLEPO, ITER
      EQUIVALENCE (W,WJ)
      DOUBLE PRECISION WJ, WK
      DATA ICALCN/0/
      DATA TYPE/'BOND       ','ANGLE      ','DIHEDRAL   '/
      DATA CALCN /'     ','ALPHA'/
      DATA NUMBRS /'0','1','2','3','4','5','6','7','8','9',' '/
      DATA FLEPO(1),FLEPO(2),FLEPO(3)/
     1' 1SCF WAS SPECIFIED, SO BFGS WAS NOT USED                 ',
     2' GRADIENTS WERE INITIALLY ACCEPTABLY SMALL                ',
     3' HERBERTS TEST WAS SATISFIED IN BFGS                      '/
      DATA FLEPO(4),FLEPO(5),FLEPO(6)/
     1' THE LINE MINIMIZATION FAILED TWICE IN A ROW.   TAKE CARE!',
     2' BFGS FAILED DUE TO COUNTS EXCEEDED. TAKE CARE!           ',
     3' PETERS TEST WAS SATISFIED IN BFGS OPTIMIZATION           '/
      DATA FLEPO(7),FLEPO(8),FLEPO(9)/
     1' THIS MESSAGE SHOULD NEVER APPEAR, CONSULT A PROGRAMMER!! ',
     2' GRADIENT TEST NOT PASSED, BUT FURTHER WORK NOT JUSTIFIED ',
     3' A FAILURE HAS OCCURRED, TREAT RESULTS WITH CAUTION!!     '/
      DATA FLEPO(10),FLEPO(11),FLEPO(12)/
     1' GEOMETRY OPTIMIZED USING NLLSQ. GRADIENT NORM MINIMIZED  ',
     2' GEOMETRY OPTIMIZED USING POWSQ. GRADIENT NORM MINIMIZED  ',
     3' CYCLES EXCEEDED, GRADIENT NOT FULLY MINIMIZED IN NLLSQ   '/
      DATA FLEPO(13),FLEPO(14),FLEPO(15)/
     1' 1SCF RUN AFTER RESTART.  GEOMETRY MIGHT NOT BE OPTIMIZED ',
     2' HEAT OF FORMATION MINIMIZED IN ONE LINE SEARCH           ',
     3' GEOMETRY OPTIMISED USING EIGENVECTOR FOLLOWING (EF).     '/
      DATA FLEPO(16)/
     1' EF-OPTIMIZED GEOMETRY.  NUMBER OF -VE ROOTS INCORRECT    '/
      DATA ITER/
     1' SCF FIELD WAS ACHIEVED                                   ',
     2'  ++++----**** FAILED TO ACHIEVE SCF. ****----++++        '/
C
C SUMMARY OF RESULTS (NOTE: THIS IS IN A SUBROUTINE SO IT
C          CAN BE USED BY THE PATH OPTION)
      SAVE
      IF(ICALCN.EQ.0)NAMFIL='**NULL**'
      IDATE=' '
      IF(IFLEPO.EQ.0) IFLEPO=7
      IUHF=MIN(INDEX(KEYWRD,' UHF'),1)+1
      PRTGRA=(INDEX(KEYWRD,' GRAD').NE.0.AND.NVAR.GT.0)
      AQCHK=(INDEX(KEYWRD,'AQUO')+INDEX(KEYWRD,'ENVAQ').NE.0)
      LINEAR=(NORBS*(NORBS+1))/2
      SINGLT=(INDEX(KEYWRD,' SING') .NE. 0)
      TRIPLT=(INDEX(KEYWRD,' TRIP') .NE. 0)
      EXCITD=(INDEX(KEYWRD,' EXCI') .NE. 0)
      CI=(INDEX(KEYWRD,' C.I.') .NE. 0)
      IF(INDEX(KEYWRD,' MINDO') .NE. 0) THEN
         CALTYP='MINDO/3'
      ELSEIF(INDEX(KEYWRD,' AM1') .NE. 0) THEN
         CALTYP='  AM1  '
      ELSEIF(INDEX(KEYWRD,' PM3') .NE. 0) THEN
         CALTYP='  PM3  '
      ELSE
         CALTYP=' MNDO  '
      ENDIF
      UHF=(IUHF.EQ.2)
      CALL DATE(IDATE)
      DEGREE=57.29577951D0
      IF(NA(1).EQ.99)THEN
         DEGREE=1.D0
         TYPE(1)='CARTESIAN X'
         TYPE(2)='CARTESIAN Y'
         TYPE(3)='CARTESIAN Z'
      ENDIF
      GNORM=0.D0
      IF(NVAR.NE.0)GNORM=SQRT(DOT(GRAD,GRAD,NVAR))
      WRITE(6,'(/,'' ----'',15(''-----''))')
      CALL WRTTXT(6)
      WRITE(6,'(//4X,A58)')FLEPO(IFLEPO)
      IITER=MAX(1,IITER)
      WRITE(6,'(4X,A58)')ITER(IITER)
      WRITE(6,'(//30X,A7,''  CALCULATION'')')CALTYP
      WRITE(6,'(55X,''VERSION '',F5.2)')VERSON
      WRITE(6,'(55X,A24)')IDATE
      IF(IITER.EQ.2)THEN
C
C   RESULTS ARE MEANINGLESS. DON'T PRINT ANYTHING!
C
         WRITE(6,'(//,'' FOR SOME REASON THE SCF CALCULATION FAILED.'',/
     1,'' THE RESULTS WOULD BE MEANINGLESS, SO WILL NOT BE PRINTED.'')')
         WRITE(6,'('' TRY TO FIND THE REASON FOR THE FAILURE BY USING ''
     1,''"PL".'',/,
     2'' CHECK YOUR GEOMETRY AND ALSO TRY USING SHIFT OR PULAY. '')')
         CALL GEOUT(1)
         STOP
      ENDIF
      IF(AQCHK) THEN
      WRITE(6,'(////10X,''FINAL HEAT OF FORMATION  '')')
      WRITE(6,'(10X,"+ DELTA-G SOLVATION     =",F13.6," KCAL")')FUNCT
      ELSE
      WRITE(6,'(////10X,''FINAL HEAT OF FORMATION ='',F13.6,'' KCAL''
     1)')FUNCT
      ENDIF
      IF(LATOM.EQ.0) WRITE(6,'(/)')
      IF(AQCHK) THEN
      WRITE(6,'(    10X,''ELECTRONIC ENERGY        '')')
      WRITE(6,'(10X,"+ DELTA-G SOLVATION     =",F13.6," EV",/)')ELECT
      ELSE
      WRITE(6,'(    10X,''ELECTRONIC ENERGY       ='',F13.6,'' EV''
     1)')ELECT
      ENDIF
      WRITE(6,'(    10X,''CORE-CORE REPULSION     ='',F17.5,'' EV''
     1)')ENUCLR
      IF(LATOM.EQ.0) WRITE(6,'(1X)')
      PRTGRA=(PRTGRA .OR. GNORM .GT. 2.D0)
      IF(PRTGRA)
     1WRITE(6,'(    10X,''GRADIENT NORM           ='',F17.5)')GNORM
      STILL=.TRUE.
      IF(LATOM.EQ.0) THEN
      IF(INDEX(KEYWRD,' AIDER').NE.0) GOTO 45
      IF(INDEX(KEYWRD,'1SCF').NE.0.AND.INDEX(KEYWRD,'GRAD').EQ.0)GOTO 45
C
C   CHECK THAT THE CARTESIAN COORDINATE GRADIENT IS ALSO SMALL
C
            IF(DOT(DXYZ,DXYZ,3*NUMAT).GT.MAX(16.D0,4*GNORM**2)
     1.AND.GNORM.LT.2.D0.AND.NCLOSE.EQ.NOPEN.AND.ID.EQ.0) THEN
               WRITE(6,'(A)')' WARNING -- GEOMETRY IS NOT AT A STATIONAR
     1Y POINT'
               STILL=.FALSE.
            ENDIF
  45  CONTINUE
      ELSE
C
C   WE NEED TO CALCULATE THE REACTION COORDINATE GRADIENT.
C
         MVAR=NVAR
         LOC11=LOC(1,1)
         LOC21=LOC(2,1)
         NVAR=1
         LOC(1,1)=LATOM
         LOC(2,1)=LPARAM
         XREACT=GEO(LPARAM,LATOM)
         CALL DERIV(GEO,GCOORD)
         NVAR=MVAR
         LOC(1,1)=LOC11
         LOC(2,1)=LOC21
         GRTYPE=' KCAL/ANGSTROM'
         IF(LPARAM.EQ.1)THEN
            WRITE(6,'(    10X,''FOR REACTION COORDINATE ='',F17.5
     1        ,'' ANGSTROMS'')')XREACT
         ELSE
            IF(NA(1).NE.99)GRTYPE=' KCAL/RADIAN  '
            WRITE(6,'(    10X,''FOR REACTION COORDINATE ='',F17.5
     1        ,'' DEGREES'')')XREACT*DEGREE
         ENDIF
         WRITE(6,'(    10X,''REACTION GRADIENT       ='',F17.5,A14
     1    )')GCOORD,GRTYPE
      ENDIF
      IF(NALPHA.GT.0)THEN
         EIONIS=-MAX(EIGS(NALPHA), EIGB(NBETA))
      ELSEIF(NELECS.EQ.1)THEN
         EIONIS=-EIGS(1)
      ELSEIF(NELECS.GT.1) THEN
         EIONIS=-MAX(EIGS(NCLOSE), EIGS(NOPEN))
      ELSE
         EIONIS=0.D0
      ENDIF
      NOPN=NOPEN-NCLOSE
C   CORRECTION TO I.P. OF DOUBLETS
      IF(NOPN.EQ.1)THEN
         I=NCLOSE*NORBS+1
         EIONIS=EIONIS+0.5D0*RJKAB(1,1)
      ENDIF
      IF(AQCHK) THEN
      WRITE(6,'(       10X,''HOMO ENERGY (EV)        ='',F13.6)')
     .-EIONIS
      ELSE
      WRITE(6,'(       10X,''IONIZATION POTENTIAL    ='',F13.6)')
     .EIONIS
      ENDIF
      WRITE(6,'(55X,A24)')IDATE
      IF( UHF ) THEN
         WRITE(6,'(      10X,''NO. OF ALPHA ELECTRONS  ='',I11)')NALPHA
         WRITE(6,'(      10X,''NO. OF BETA  ELECTRONS  ='',I11)')NBETA
      ELSE
         WRITE(6,'(      10X,''NO. OF FILLED LEVELS    ='',I11)')NCLOSE
         IF(NOPN.NE.0) THEN
            WRITE(6,'(   10X,''AND NO. OF OPEN LEVELS  ='',I11)')NOPN
         ENDIF
      ENDIF
      SUMW=0
      DO 10 I=1,NUMAT
   10 SUMW=SUMW+ATMASS(I)
      IF(SUMW.GT.0.1D0)
     1WRITE(6,'(    10X,''MOLECULAR WEIGHT        ='',F11.3)')SUMW
      IF(LATOM.EQ.0) WRITE(6,'(/)')
      WRITE(6,'(10X,''SCF CALCULATIONS  =   '',I14 )') NSCF
      TIM=SECOND()-TIME0
      I=TIM*0.000001D0
      TIM=TIM-I*1000000
      CALL TIMOUT(6,TIM)
      IF( NDEP .NE. 0 )CALL SYMTRY
      DO 20 I=1,NVAR
   20 XPARAM(I)=GEO(LOC(2,I),LOC(1,I))
      CALL GMETRY(GEO,COORD)
      IF(PRTGRA)THEN
         WRITE(6,'(///7X,''FINAL  POINT  AND  DERIVATIVES'',/)')
         WRITE(6,'(''   PARAMETER     ATOM    TYPE  ''
     1    ,''          VALUE       GRADIENT'')')
      ENDIF
      SUM=0.5D0
      DO 30 I=1,NUMAT
   30 SUM=SUM+CORE(NAT(I))
      I=SUM
      KCHRGE=I-NCLOSE-NOPEN-NALPHA-NBETA
C
C    WRITE OUT THE GEOMETRIC VARIABLES
C
      IF(PRTGRA) THEN
         DO 40 I=1,NVAR
            J=LOC(2,I)
            K=LOC(1,I)
            L=LABELS(K)
            XI=XPARAM(I)
            IF(J.NE.1) XI=XI*DEGREE
            IF(J.EQ.1.OR.NA(1).EQ.99)THEN
               GTYPE='KCAL/ANGSTROM'
            ELSE
               GTYPE='KCAL/RADIAN  '
            ENDIF
   40    WRITE(6,'(I7,I11,1X,A2,4X,A11,F13.6,F13.6,2X,A13)')
     1I,K,ELEMNT(L),TYPE(J),XI,GRAD(I),GTYPE
      ENDIF
C
C     WRITE OUT THE GEOMETRY
C
      WRITE(6,'(///)')
      CALL GEOUT(1)
      IF (INDEX(KEYWRD,' NOINTER') .EQ. 0) THEN
C
C   WRITE OUT THE INTERATOMIC DISTANCES
C
         L=0
         DO 50 I=1,NUMAT
            DO 50 J=1,I
               L=L+1
   50    RXYZ(L)=SQRT((COORD(1,I)-COORD(1,J))**2+
     1                         (COORD(2,I)-COORD(2,J))**2+
     2                         (COORD(3,I)-COORD(3,J))**2)
         WRITE(6,'(//10X,''  INTERATOMIC DISTANCES'')')
         CALL VECPRT(RXYZ,NUMAT)
      ENDIF
      DO 60 I=1,NORBS
   60 IF(EIGS(I).LT.-999.D0.OR.EIGS(I).GT.1000.D0)EIGS(I)=0.D0
      DO 70 I=1,NORBS
   70 IF(EIGB(I).LT.-999.D0.OR.EIGB(I).GT.1000.D0)EIGS(I)=0.D0
      IF(ISYBYL.EQ.1) THEN
C
C  THE FOLLOWING OPEN STATEMENTS ARE NON-STANDARD.  IF THIS CAUSES 
C  DIFFICULTY REPLACE THEM WITH
      OPEN(UNIT=16,FILE=GETNAM('FOR016'),STATUS='NEW',ERR=31)
      GOTO 32
  31  OPEN(UNIT=16,FILE=GETNAM('FOR016'),STATUS='OLD')
      WRITE(6,'(A)') 'Error opening SYBYL MOPAC output'
  32  CONTINUE
C#      OPEN(UNIT=16,FILE=GETNAM('FOR016'),CARRIAGECONTROL='LIST',
C#     +STATUS='NEW',ERR=31)
C#      GOTO 32
C#  31  OPEN(UNIT=16,FILE=GETNAM('FOR016'),CARRIAGECONTROL='LIST',
C#     +STATUS='OLD')
C#      WRITE(6,'(A)') 'Error opening SYBYL MOPAC output'
C#  32  CONTINUE
      ENDIF
      IF(NORBS.GT.0)THEN
         IF (INDEX(KEYWRD,' VECT') .NE. 0) THEN
            WRITE(6,'(//10X,A5,'' EIGENVECTORS  '')')CALCN(IUHF)
            CALL MATOUT (C,EIGS,NORBS,NORBS,NORBS)
            IF(UHF) THEN
               WRITE(6,'(//10X,'' BETA EIGENVECTORS  '')')
               CALL MATOUT (CBETA,EIGB,NORBS,NORBS,NORBS)
            ENDIF
         ELSE
            WRITE(6,'(//10X,A5,''   EIGENVALUES'',/)')CALCN(IUHF)
            WRITE(6,'(8F10.5)')(EIGS(I),I=1,NORBS)
            IF(UHF) THEN
               WRITE(6,'(//10X,'' BETA EIGENVALUES '')')
               WRITE(6,'(8F10.5)')(EIGB(I),I=1,NORBS)
            ENDIF
         ENDIF
      ENDIF
      WRITE(6,'(//13X,'' NET ATOMIC CHARGES AND DIPOLE '',
     1''CONTRIBUTIONS'',/)')
      WRITE(6,'(8X,'' ATOM NO.   TYPE          CHARGE        ATOM''
     1,''  ELECTRON DENSITY'')')
      CALL CHRGE(P,Q)
      DO 80 I=1,NUMAT
         L=NAT(I)
         Q2(I)=CORE(L) - Q(I)
   80 WRITE(6,'(I12,9X,A2,4X,F13.4,F16.4)')
     1I,ELEMNT(L),Q2(I),Q(I)
      DIP= DIPOLE(P,Q2,COORD,DUMY,1)
      IF (INDEX(KEYWRD,' NOXYZ') .EQ. 0) THEN
         WRITE(6,'(//10X,''CARTESIAN COORDINATES '',/)')
         WRITE(6,'(4X,''NO.'',7X,''ATOM'',15X,''X'',
     1  9X,''Y'',9X,''Z'',/)')
         WRITE(6,'(I6,8X,A2,14X,3F10.4)')
     1  (I,ELEMNT(NAT(I)),(COORD(J,I),J=1,3),I=1,NUMAT)
      ENDIF
      IF(NORBS.GT.0) THEN
         IF (INDEX(KEYWRD,' K=') .NE. 0)THEN
C
C  GO INTO BRILLOUIN ZONE MODE
C
            I=INDEX(KEYWRD,' K=')
            STEP=READA(KEYWRD,I)
            MONO3=NLAST(NINT(READA(KEYWRD(I:),INDEX(KEYWRD(I:),','))))
         IF(UHF)WRITE(6,'(A)')'  ALPHA BANDS'
            CALL BRLZON(F, FMAT2D, NORBS, SEC, VEC, ALBAND,MONO3,STEP,2)
         IF(UHF)THEN
         WRITE(6,'(A)')'  BETA BANDS'
         CALL BRLZON(FB, FMAT2D, NORBS, SEC, VEC, ALBAND,MONO3,STEP,2)
         ENDIF
         ENDIF
         IF(ISYBYL.EQ.1)THEN
            NFILLD=MAX(NCLOSE,NALPHA,NBETA)
            CALL MPCSYB(NUMAT,COORD,Q2,1,EIGS,NFILLD,FUNCT,EIONIS
     1                 ,KCHRGE,DIP)
         ENDIF
         IF (INDEX(KEYWRD,' FOCK') .NE. 0) THEN
            WRITE(6,'('' FOCK MATRIX IS '')')
            CALL VECPRT(F,NORBS)
         ENDIF
         IF (INDEX(KEYWRD,' DENS') .NE. 0) THEN
            WRITE(6,'(//,20X,'' DENSITY MATRIX IS '')')
            CALL VECPRT(P,NORBS)
         ELSE
            WRITE(6,'(//10X,''ATOMIC ORBITAL ELECTRON POPULATIONS'',/)')
            WRITE(6,'(8F10.5)')(P((I*(I+1))/2),I=1,NORBS)
         ENDIF
         IF(INDEX(KEYWRD,' PI') .NE. 0) THEN
            WRITE(6,'(//10X,''SIGMA-PI BOND-ORDER MATRIX'')')
            CALL DENROT
         ENDIF
         IF(UHF) THEN
            SZ=ABS(NALPHA-NBETA)*0.5D0
            SS2=SZ*SZ
            L=0
            DO 100 I=1,NORBS
               DO 90 J=1,I
                  L=L+1
                  PA(L)=PA(L)-PB(L)
   90          SS2=SS2+PA(L)**2
  100       SS2=SS2-0.5D0*PA(L)**2
            WRITE(6,'(//20X,''(SZ)    ='',F10.6)')SZ
            WRITE(6,'(  20X,''(S**2)  ='',F10.6)')SS2
            IF(INDEX(KEYWRD,' SPIN') .NE. 0) THEN
               WRITE(6,'(//10X,''SPIN DENSITY MATRIX'')')
               CALL VECPRT(PA,NORBS)
            ELSE
               WRITE(6,'(//10X,''ATOMIC ORBITAL SPIN POPULATIONS'',/)')
               WRITE(6,'(8F10.5)')(PA((I*(I+1))/2),I=1,NORBS)
            ENDIF
            IF(INDEX(KEYWRD,' HYPERFINE') .NE. 0) THEN
C
C  WORK OUT THE HYPERFINE COUPLING CONSTANTS.
C
               WRITE(6,'(//10X,''    HYPERFINE COUPLING COEFFICIENTS'',/
     1)')
               J=(NALPHA-1)*NORBS
               DO 110 K=1,NUMAT
                  I=NFIRST(K)
C#          WRITE(6,'('' PA:'',F13.6,'' C('',I2,''+'',I3,''):'',
C#     +F13.5)')PA((I*(I+1))/2),I,J,C(I+J)
  110          Q(K)=PA((I*(I+1))/2)*0.3333333D0+C(I+J)**2*0.66666666D0
               WRITE(6,'(5(2X,A2,I2,F9.5,1X))')
     1    (ELEMNT(NAT(I)),I,Q(I),I=1,NUMAT)
            ENDIF
            DO 120 I=1,LINEAR
  120       PA(I)=P(I)-PB(I)
         ENDIF
         IF (INDEX(KEYWRD,' BONDS') .NE. 0) THEN
            IF(NBETA.EQ.0)THEN
               WRITE(6,'(/10X,''BONDING CONTRIBUTION OF EACH M.O.'',/)')
               CALL MOLVAL(C,P,NORBS,2.D0)
            ELSE
               WRITE(6,'(/10X,''BONDING CONTRIBUTION OF EACH ALPHA M.O.'
     1',/)')
               CALL MOLVAL(C,P,NORBS,1.D0)
               WRITE(6,'(/10X,''BONDING CONTRIBUTION OF EACH BETA  M.O.'
     1',/)')
               CALL MOLVAL(C,P,NORBS,1.D0)
            ENDIF
            CALL BONDS(P)
         ENDIF
         I=NCLOSE+NALPHA
         IF (INDEX(KEYWRD,' LOCAL') .NE. 0) THEN
            CALL LOCAL(C,NORBS,I,EIGS)
            IF(NBETA.NE.0)THEN
               WRITE(6,'(//10X,'' LOCALIZED BETA MOLECULAR ORBITALS'')')
               CALL LOCAL(CBETA,NORBS,NBETA,EIGB)
            ENDIF
         ENDIF
         IF (INDEX(KEYWRD,' 1ELE') .NE. 0) THEN
            WRITE(6,'('' FINAL ONE-ELECTRON MATRIX '')')
            CALL VECPRT(H,NORBS)
         ENDIF
         IF(INDEX(KEYWRD,' ENPART') .NE. 0)
     1CALL ENPART(UHF,H,PA,PB,P,Q,COORD)
      ENDIF
      DO 130 I=1,107
  130 NELEMT(I)=0
      DO 140 I=1,NUMAT
         IGO=NAT(I)
         IF (IGO.GT.107) GO TO 140
         NELEMT(IGO)=NELEMT(IGO)+1
  140 CONTINUE
      ICHFOR=0
      IF (NELEMT(6).EQ.0) GO TO 150
      ICHFOR=1
      IELEMT(1)=ELEMNT(6)
      NZS=NELEMT(6)
      IF (NZS.LT.10) THEN
         IF (NZS.EQ.1) THEN
            IEL1(1)=11
         ELSE
            IEL1(1)=NZS+1
         ENDIF
         IEL2(1)=11
      ELSE
         KFRST=NZS/10
         KSEC=NZS-(10*KFRST)
         IEL1(1)=KFRST+1
         IEL2(1)=KSEC+1
      ENDIF
  150 NELEMT(6)=0
      DO 160 I=1,107
         IF (NELEMT(I).EQ.0) GO TO 160
         ICHFOR=ICHFOR+1
         IELEMT(ICHFOR)=ELEMNT(I)
         NZS=NELEMT(I)
         IF (NZS.LT.10) THEN
            IF (NZS.EQ.1) THEN
               IEL1(ICHFOR)=11
            ELSE
               IEL1(ICHFOR)=NZS+1
            ENDIF
            IEL2(ICHFOR)=11
         ELSE
            KFRST=NZS/10
            KSEC=NZS-(10*KFRST)
            IEL1(ICHFOR)=KFRST+1
            IEL2(ICHFOR)=KSEC+1
         ENDIF
  160 CONTINUE
      IF(INDEX(KEYWRD,' DENOUT') .NE. 0) THEN
         OPEN(UNIT=10,FILE=GETNAM('FOR010'),
     +STATUS='UNKNOWN',FORM='UNFORMATTED')
         REWIND 10
         WRITE(10)(PA(I),I=1,LINEAR)
         IF(UHF)WRITE(10)(PB(I),I=1,LINEAR)
         CLOSE (10)
      ENDIF
      IF((CI.OR.NOPEN.NE.NCLOSE.AND.FRACT.NE.2.D0.AND.FRACT.NE.0.D0
     1 .OR.INDEX(KEYWRD,' SIZE').NE.0)
     2 .AND. INDEX(KEYWRD,' MECI')+INDEX(KEYWRD,' ESR').NE.0)THEN
         WRITE(6,'(//10X,
     1''MULTI-ELECTRON CONFIGURATION INTERACTION CALCULATION'',//)')
         LAST=3
         X=MECI(EIGS,C)
      ENDIF
      IF (INDEX(KEYWRD,' MULLIK') +INDEX(KEYWRD,' GRAPH') .NE. 0) THEN
         IF (INDEX(KEYWRD,' MULLIK') .NE. 0)
     1   WRITE(6,'(/10X,'' MULLIKEN POPULATION ANALYSIS'')')
         CALL MULLIK(C,H,F,NORBS,P,RXYZ)
         IF (INDEX(KEYWRD,' GRAPH') .NE. 0)
     1   WRITE(6,'(/10X,'' DATA FOR GRAPH WRITTEN TO DISK'')')
      ENDIF
      IF((INDEX(KEYWRD,'AQUO')+INDEX(KEYWRD,'ENVAQ')).EQ.0) GO TO 2156
      WRITE(6,'(//,"f(GB) factors as defined by Still et al., J. Amer. C
     1hem. Soc. 1990, 112, 6127")')
      CALL VECPRT(FGB,NUMAT)
      WRITE(6,'(//)')
      WRITE(6,'(" Generalized Born Polarization Energy decomposition"/)'
     1)
      DO 260 I=1,NUMAT
      CELEID=1.D0/78.3D0
      AK=-166.D0*(1.D0-CELEID)
      JP=I+1
      BP(I)=0.D0
      DO 270 J=1,I
      IJ=(I*(I-1))/2+J
      EBP(IJ)=AK*QEFF(I)*QEFF(J)/FGB(IJ)
270   BP(I)=BP(I)+EBP(IJ)
      DO 280 K=JP,NUMAT
      IJ=(K*(K-1))/2+I
      EBP(IJ)=AK*QEFF(I)*QEFF(K)/FGB(IJ)
280   BP(I)=BP(I)+EBP(IJ)
260   CONTINUE
      CALL VECPRT(EBP,NUMAT)
      WRITE(6,'(//" BY ATOM:"/)')
      DO 290 L=1,NUMAT
290   WRITE(6,'(" ATOM #",I3," BORN ENERGY ",F7.2," KCAL/MOL")')L,BP(L)
      WRITE(6,'(/)')
      DO 2151 I=1,NUMAT
2151  WRITE(6,2152)I,ATAR(I),SURFACT(I)
2152  FORMAT(' FOR ATOM #',I3,' WITH AREA',F7.2,' ANG^2 THE SURFACE CORR
     1ECTION IS',F7.2,' KCAL/MOL')
      WRITE(6,'(//)')
      WRITE(6,'(" BY ELEMENT",/)')
      BPGT=0.D0
      DO 300 M=1,100
      BPT=0.D0
      ST=0.D0
      EST(M)=0.D0
      IF(INDEX(KEYWRD,'ENVAQ').EQ.0) GO TO 305
      DO 307 N=1,NUMAT
      IF(ITYPE(N).NE.M) GO TO 307
      EST(M)=EST(M)+ATAR(N)
307   CONTINUE
305   DO 310 N=1,NUMAT
      IF(NAT(N).NE.M) GO TO 310
      BPT=BPT+BP(N)
      ST=ST+SURFACT(N)
      BPGT=BPGT+BP(N)
310   CONTINUE
      IF(ST.EQ.0.D0.AND.BPT.EQ.0.D0) GO TO 300
      WRITE(6,'(" ATOMIC #",I3," BORN: ",F7.2," SURFACE: ",F7.2," TOTAL:
     1 ",F7.2," KCAL/MOL")')M,BPT,ST,BPT+ST
300   CONTINUE
      WRITE(6,'(/," ** NET ENERGIES:  ",F7.2,10X,F7.2,8X,F7.2," KCAL/MOL
     1")')BPGT,SURFCT,BPGT+SURFCT
      WRITE(6,'(//,"**** NOTA BENE:  THE NET SOLVATION ENERGY IS FOR THI
     1S EXACT MOLECULAR",/,"GEOMETRY! The true AM1-SM1 or AM1-SM1a solva
     2tion energy",/,"should be obtained as the difference between the h
     3eat of",/,"formation for the relaxed aqueous system and that for t
     4he",/,"relaxed gas-phase system.",/)')
      IF(INDEX(KEYWRD,'NOFUL').NE.0) WRITE(6,306)
306   FORMAT(' * THIS IS NOT A TRUE STATIONARY POINT ')
      IF(INDEX(KEYWRD,'ENVAQ').EQ.0) GO TO 312
      DO 311 N=1,18
      IF(EST(N).EQ.0.D0) GO TO 311
      WRITE(6,'(" TOTAL AREA FOR ENVIRONMENT TYPE ",I2," IS ",F7.2," A^2
     1")')N,EST(N)
 311   CONTINUE
312   IF(INDEX(KEYWRD,'FOCK').EQ.0) GO TO 2156
      WRITE(6,'(///)')
      DO 2154 J=1,NUMAT
2154  WRITE(6,2155)J,BP(J)
2155  FORMAT('ATOM #',I3,' BORN CONTRIBUTION TO FOCK DIAGONAL:',F12.6,
     1' EV')
C
C   *****************************************************************
C   *                                                               *
C   *      SUMMARY OF OUTPUT ON FILE IWRITE ( . ARC FILE )          *
C   *                                                               *
C   *****************************************************************
C
C     NOT DONE IF OPTIMISATION NOT ACHIEVED
2156  IF (IFLEPO.EQ.5 .OR. IFLEPO.EQ.9 .OR. IFLEPO.EQ.12) RETURN
C
C  NOTE THAT THE DENSITY, H AND F MATRICES ARE CORRUPTED BY A
C  CALL TO MULLIK.
      IF(ISYBYL.EQ.1) THEN
         IF (INDEX(KEYWRD,'MULLIK').EQ.0) THEN
            CALL MPCPOP(C,0)
         ELSE
            CALL MPCPOP(C,1)
         ENDIF
         CLOSE(16)
      ENDIF
      IF(ICALCN.NE.NUMCAL)THEN
         IF(NAMFIL.EQ.'**NULL**') THEN
         NAMFIL=GETNAM('FOR012')
         INAM=ICHAR('a')
         JNAM=INAM
         JEND=INDEX(NAMFIL,' ')
         IEND=JEND+1
         ENDIF
  162    CLOSE (12)
         OPEN(UNIT=12,FILE=NAMFIL,STATUS='NEW',ERR=163)
         GOTO 164
  163    NAMFIL(IEND:IEND)=CHAR(INAM)
         NAMFIL(JEND:JEND)=CHAR(JNAM)
         IF(INAM.EQ.ICHAR('z'))THEN
         INAM=INAM-26
         JNAM=JNAM+1
         ENDIF
         INAM=INAM+1
         GOTO 162
  164    REWIND 12
         ICALCN=NUMCAL
      ENDIF
      IWRITE=12
  170 WRITE(IWRITE,'(//20X,'' SUMMARY OF '',A7,
     1'' CALCULATION'',/)')CALTYP
      WRITE(IWRITE,'(60X,''VERSION '',F5.2)')VERSON
      WRITE (IWRITE,180) (IELEMT(I),NUMBRS(IEL1(I)),NUMBRS(IEL2(I))
     1,I=1,ICHFOR)
  180 FORMAT (//,1X,17(A2,A1,A1))
      WRITE(IWRITE,'(55X,A24)')IDATE
      CALL WRTTXT(IWRITE)
      WRITE(IWRITE,'(//4X,A58)')FLEPO(IFLEPO)
      WRITE(IWRITE,'(4X,A58)')ITER(IITER)
      IXRT=IWRITE
      IF(AQCHK) THEN
      WRITE(IXRT,'(/10X,''FINAL HEAT OF FORMATION  '')')
      WRITE(IXRT,'(10X,"+ DELTA-G SOLVATION     =",F13.6," KCAL")')FUNCT
      ELSE
      WRITE(IXRT,'(/10X,''FINAL HEAT OF FORMATION ='',F13.6,'' KCAL''
     1)')FUNCT
      ENDIF
      IF(AQCHK) THEN
      WRITE(IXRT,'(    10X,''ELECTRONIC ENERGY        '')')
      WRITE(IXRT,'(10X,"+ DELTA-G SOLVATION     =",F13.6," EV")')ELECT
      ELSE
      WRITE(IXRT,'(    10X,''ELECTRONIC ENERGY       ='',F13.6,'' EV''
     1)')ELECT
      ENDIF
      WRITE(IWRITE,'(  10X,''CORE-CORE REPULSION     =''
     1,F17.6,'' EV'')')ENUCLR
      IF(PRTGRA)
     1WRITE(IWRITE,'(  10X,''GRADIENT NORM           =''
     2,F17.6)')GNORM
      IF(LATOM.EQ.0) THEN
         IF(.NOT.STILL) WRITE(IWRITE,'(A)')
     1' WARNING -- GEOMETRY IS NOT AT A STATIONARY POINT'
      ELSE
         GRTYPE=' KCAL/ANGSTROM'
         IF(LPARAM.EQ.1)THEN
            WRITE(IWRITE,'(    10X,''FOR REACTION COORDINATE ='',F17.4
     1        ,'' ANGSTROMS'')')XREACT
         ELSE
            IF(NA(1).NE.99)GRTYPE=' KCAL/RADIAN  '
            WRITE(IWRITE,'(    10X,''FOR REACTION COORDINATE ='',F17.4
     1        ,'' DEGREES'')')XREACT*DEGREE
         ENDIF
         WRITE(IWRITE,'(    10X,''REACTION GRADIENT       ='',F17.6,A14
     1    )')GCOORD,GRTYPE
      ENDIF
      WRITE(IWRITE,'(  10X,''DIPOLE                  =''
     1,F16.5, '' DEBYE'')')DIP
      IF(UHF) THEN
         WRITE(IWRITE,'(  10X,''(SZ)                    ='',F17.6)')SZ
         WRITE(IWRITE,'(  10X,''(S**2)                  ='',F17.6)')SS2
         WRITE(IWRITE,'(  10X,''NO. OF ALPHA ELECTRONS  ='',I10)')NALPHA
         WRITE(IWRITE,'(  10X,''NO. OF BETA  ELECTRONS  ='',I10)')NBETA
      ELSE
         WRITE(IWRITE,'(  10X,''NO. OF FILLED LEVELS    ='',I10)')NCLOSE
         NOPN=NOPEN-NCLOSE
         IF(NOPN.NE.0)
     1WRITE(IWRITE,'(  10X,''AND NO. OF OPEN LEVELS  ='',I10)')NOPN
      ENDIF
      IF(CI)
     1WRITE(IWRITE,'(  10X,''CONFIGURATION INTERACTION WAS USED'')')
      IF(KCHRGE.NE.0)
     1WRITE(IWRITE,'(  10X,''CHARGE ON SYSTEM        ='',I10)')KCHRGE
      IF(AQCHK) THEN
      WRITE(IWRITE,'(       10X,''HOMO ENERGY (EV)        ='',F13.6)')
     .-EIONIS
      ELSE
      WRITE(IWRITE,'(       10X,''IONIZATION POTENTIAL    ='',F13.6)')
     .EIONIS
      ENDIF
      WRITE(IWRITE,'(  10X,''MOLECULAR WEIGHT        ='',F14.3)')SUMW
      WRITE(IWRITE,'(  10X,''SCF CALCULATIONS        =''
     1,I10)') NSCF
      TIM=SECOND()-TIME0
      CALL TIMOUT(IWRITE,TIM)
      WRITE(IWRITE,'(//10X,''FINAL GEOMETRY OBTAINED'',36X,''CHARGE'')')
      CALL GEOUT(IWRITE)
      IF(INDEX(KEYWRD,' AIGOUT').NE.0)THEN
         WRITE(IWRITE,'(//,A)')'  GEOMETRY IN GAUSSIAN Z-MATRIX STYLE'
         CALL WRTTXT(IWRITE)
         CALL GEOUTG(IWRITE)
      ENDIF
      IF(IWRITE.NE.11.AND.INDEX(KEYWRD,' NOLOG').EQ.0)THEN
         IWRITE=11
         GOTO 170
      ENDIF
      NSCF=0
      RETURN
      END
      SUBROUTINE TIMOUT(NOUT,TIM)
C
C     CONVERT THE TIME FROM SECONDS TO DAYS, HOURS, MINUTES, AND SECONDS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION MINS, MINPHR
C
C
      DATA HRSPD /24.0D0/,    MINPHR /60.0D0/
      DATA SECPD /86400.0D0/, SECPMI /60.0D0/
C
      DAYS = TIM / SECPD
      IDAYS = INT(DAYS)
      HOURS = (DAYS - FLOAT(IDAYS)) * HRSPD
      IHOURS = INT(HOURS)
      MINS = (HOURS - FLOAT(IHOURS)) * MINPHR
      IMINS = INT(MINS)
      SECS = (MINS - FLOAT(IMINS)) * SECPMI
C
      IF (IDAYS .GT. 1) THEN
         WRITE (NOUT,10) IDAYS,IHOURS,IMINS,SECS
      ELSE IF (IDAYS .EQ. 1) THEN
         WRITE (NOUT,20) IDAYS,IHOURS,IMINS,SECS
      ELSE IF (IHOURS .GT. 0) THEN
         WRITE (NOUT,30) IHOURS,IMINS,SECS
      ELSE IF (IMINS .GT. 0) THEN
         WRITE (NOUT,40) IMINS,SECS
      ELSE
         WRITE (NOUT,50) SECS
      END IF
C
   10 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'DAYS',2X,I2,1X,'HOURS',
     1        1X,I2,1X,'MINUTES AND',1X,F7.3,1X,'SECONDS')
   20 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'DAY',2X,I2,1X,'HOURS',
     1        1X,I2,1X,'MINUTES AND',1X,F7.3,1X,'SECONDS')
   30 FORMAT (10X,'COMPUTATION TIME = 'I2,1X,'HOURS',
     1        1X,I2,1X,'MINUTES AND',1X,F7.3,1X,'SECONDS')
   40 FORMAT (10X,'COMPUTATION TIME = ',I2,1X,'MINUTES AND',
     1        1X,F7.3,1X,'SECONDS')
   50 FORMAT (10X,'COMPUTATION TIME = ',F7.3,1X,'SECONDS')
      END
      SUBROUTINE MPCPOP(C,ICOK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
C
C This subroutine calculates the total Mulliken populations on the
C   atoms by summing the diagonal elements from the  Mulliken
C   population analysis.
C
      COMMON / MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /CORE/ CORE(107)
      DIMENSION C(MORB2),POP(NUMATM),CHRG(NUMATM)
      WRITE(16,'(I4,5X'' MULLIKEN POPULATION AND CHARGE'')',ERR=40)ICOK
C
C ICOK = 1 ==> PRINT POPULATIONS
C ICOK = 0 ==> KEYWORD mulliken = .f.
C         NO POPULATION ANALYSIS PERFORMED
C
      IF (ICOK.NE.0) THEN
         DO 20 I = 1,NUMAT
            IF = NFIRST(I)
            IL = NLAST(I)
            SUM = 0.0
            POP(I) = 0.0
            CHRG(I) = 0.0
            DO 10 J = IF,IL
C
C    Diagonal element of mulliken matrix
C
               SUM = SUM + C((J*(J+1))/2)
   10       CONTINUE
            K = NAT(I)
C
C    Mulliken population for i'th atom
C
            POP(I) = SUM
            CHRG(I) = CORE(K) - POP(I)
   20    CONTINUE
         WRITE(6,'(///10X,''MULLIKEN POPULATIONS AND CHARGES'')')
         DO 30 J = 1,NUMAT
            WRITE(6,60) J, POP(J), CHRG(J)
            WRITE(16,70,ERR=40) POP(J), CHRG(J)
   30    CONTINUE
      ENDIF
      RETURN
   40 WRITE(6,'(A)') 'Error writing SYBYL Mulliken population output'
      RETURN
   50 FORMAT(//,5X,'ATOM',8X,'POPULATION',6X,'CHARGE')
   60 FORMAT(5X,I4,4X,F11.6,6X,F11.6)
   70 FORMAT(2F12.6)
      END
C
C This subroutine writes out the optimized geometry and atomic charges
C   for a MOPAC run.
C
      SUBROUTINE MPCSYB(NUMAT,COORD,CHR,ICOK,EIGS,NCLOSE,FUNCT
     1                       ,EIONIS,KCHRGE,DIP)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION COORD(3, NUMAT), CHR(NUMAT),EIGS(MAXORB)
C  Write out the charge flag and number of atoms
      WRITE(16,'(2I4)', ERR=30) ICOK,NUMAT
C  Write out the coordinates and charges
      DO 10 I=1, NUMAT
         WRITE(16,'(4F12.6)', ERR=30) (COORD(J, I), J=1, 3), CHR(I)
   10 CONTINUE
      I1 = MAX(1,NCLOSE - 1)
      I2 = MIN(MAXORB,NCLOSE + 2)
C
C  Write out the 2 highest and 2 lowest orbital energies
C
      WRITE(16,20,ERR=30)(EIGS(J),J=I1,I2),NCLOSE
   20 FORMAT(4F12.6,2X,I4,2X,'HOMOs,LUMOs,# of occupied MOs')
C
C  Write out the Heat of Formation and Ionisation Potential
C
      WRITE(16,'(2F12.6,4X,''HF and IP'')',ERR=30) FUNCT,EIONIS
C
C  Write out the Dipole Moment
C
      IF(KCHRGE.NE.0) DIP = 0.0
      WRITE(16, '(I4,F10.3,''  Charge,Dipole Moment'')', ERR=30)
     1KCHRGE, DIP
      RETURN
   30 WRITE(6,'(A)') 'Error writing SYBYL MOPAC output'
      RETURN
      END
