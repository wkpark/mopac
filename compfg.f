      SUBROUTINE COMPFG(XPARAM,INT,ESCF,FULSCF,GRAD,LGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION XPARAM(MAXPAR),GRAD(MAXPAR)
      LOGICAL LGRAD, FULSCF
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR),IDUMY,DUMY(MAXPAR)
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /ATHEAT/ ATHEAT
      COMMON /WMATRX/ WJ(N2ELEC), WK(N2ELEC)
      COMMON /ENUCLR/ ENUCLR
      COMMON /NATYPE/ NZTYPE(107),MTYPE(30),LTYPE
      COMMON /ELECT / ELECT
      PARAMETER (MDUMY=MAXPAR**2-MPACK)
      COMMON /SCRACH/ RXYZ(MPACK), XDUMY(MDUMY)
      COMMON /HMATRX/ H(MPACK)
      COMMON /MOLMEC/ HTYPE(4),NHCO(4,20),NNHCO,ITYPE,USEMM
     1       /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     2                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     3                NCLOSE,NOPEN,NDUMY,FRACT
C***********************************************************************
C
C   COMPFG CALCULATES (A) THE HEAT OF FORMATION OF THE SYSTEM, AND
C                     (B) THE GRADIENTS, IF LGRAD IS .TRUE.
C
C   ON INPUT  XPARAM = ARRAY OF PARAMETERS TO BE USED IN INTERNAL COORDS
C             LGRAD  = .TRUE. IF GRADIENTS ARE NEEDED, .FALSE. OTHERWISE
C             INT    = NEVER USED, RESERVED FOR FUTURE USE
C             FULSCF = .TRUE. IF FULL SCF TO BE DONE, .FALSE. OTHERWISE.
C
C   ON OUTPUT ESCF  = HEAT OF FORMATION.
C             GRAD   = ARRAY OF GRADIENTS, IF LGRAD = .TRUE.
C
C***********************************************************************
      COMMON /KEYWRD/KEYWRD
      CHARACTER*80 KEYWRD
      REAL WJ, WK
      LOGICAL FIRST, DEBUG, INT, PRINT, ANALYT, LARGE, USEMM
      DIMENSION COORD(3,NUMATM), W(N2ELEC), DEGREE(3)
      EQUIVALENCE (W,WJ)
      DATA FIRST /.TRUE./
C                 MNDO     AM1      PM3      MINDO/3
      DATA HTYPE/ 6.1737D0,  3.3191D0,  7.1853D0,  1.7712D0/
      IF (FIRST) THEN
         FIRST=.FALSE.
C#         WRITE(6,'(A,F12.4)')' MM CORRECTION:',HTYPE(ITYPE)
         LTYPE=0
         DO 30 I=1,NUMAT
            IF(NAT(I).LT.99)THEN
               DO 10 J=1,LTYPE
   10          IF(NAT(I).EQ.MTYPE(J)) GOTO 20
               LTYPE=LTYPE+1
               MTYPE(LTYPE)=NAT(I)
               NZTYPE(NAT(I))=LTYPE
C
C       LTYPE = NUMBER OF TYPES OF REAL ATOM PRESENT
C       MTYPE = TYPES OF REAL ATOMS PRESENT
               J=LTYPE
   20          CONTINUE
            ENDIF
   30    CONTINUE
         ANALYT=(INDEX(KEYWRD,'ANALYT').NE.0)
         IF(ANALYT)CALL SETUPG
         DEGREE(1)=1.D0
         IF(INDEX(KEYWRD,' XYZ').NE.0)THEN
            DEGREE(2)=1.D0
         ELSE
            DEGREE(2)=180.D0/3.141592652589D0
         ENDIF
         DEGREE(3)=DEGREE(2)
         LARGE=(INDEX(KEYWRD,'LARGE') .NE. 0)
         PRINT=(INDEX(KEYWRD,'COMPFG') .NE. 0)
         DEBUG=(INDEX(KEYWRD,'DEBUG') .NE. 0 .AND. PRINT)
      ENDIF
C
C SET UP COORDINATES FOR CURRENT CALCULATION
C
C       PLACE THE NEW VALUES OF THE VARIABLES IN THE ARRAY GEO.
C       MAKE CHANGES IN THE GEOMETRY.
      DO 40 I=1,NVAR
         K=LOC(1,I)
         L=LOC(2,I)
   40 GEO(L,K)=XPARAM(I)
C#      WRITE(6,'(3F18.11)')(XPARAM(I),I=1,3)
C      IMPOSE THE SYMMETRY CONDITIONS + COMPUTE THE DEPENDENT-PARAMETERS
      IF(NDEP.NE.0) CALL SYMTRY
C      NOW COMPUTE THE ATOMIC COORDINATES.
      IF( DEBUG ) THEN
         IF( LARGE ) THEN
            K=NUMAT
         ELSE
            K=MIN(5,NUMAT)
         ENDIF
         WRITE(6,FMT='('' INTERNAL COORDS'',/100(/,3F12.6))')
     1            ((GEO(J,I)*DEGREE(J),J=1,3),I=1,K)
      END IF
      CALL GMETRY(GEO,COORD)
      IF( DEBUG ) THEN
         IF( LARGE ) THEN
            K=NUMAT
         ELSE
            K=MIN(5,NUMAT)
         ENDIF
         WRITE(6,FMT='('' CARTESIAN COORDS'',/100(/,3F12.6))')
     1            ((COORD(J,I),J=1,3),I=1,K)
      END IF
      IF(ANALYT)REWIND 2
      CALL HCORE(COORD, H, W, WJ, WK, ENUCLR)
C
C COMPUTE THE HEAT OF FORMATION.
C
      IF(NORBS.GT.0.AND.NELECS.GT.0) THEN
         CALL ITER(H, W, WJ, WK, ELECT, FULSCF,.TRUE.)
      ELSE
         ELECT=0.D0
      ENDIF
      ESCF=(ELECT+ENUCLR)*23.061D0+ATHEAT
      DO 50 I=1,NNHCO
         CALL DIHED(COORD,NHCO(1,I),NHCO(2,I),NHCO(3,I),NHCO(4,I),ANGLE)
         ESCF=ESCF+HTYPE(ITYPE)*SIN(ANGLE)**2
   50 CONTINUE
      IF(PRINT) WRITE(6,'(/10X,'' HEAT OF FORMATION'',G30.17)')ESCF
C
C FIND DERIVATIVES IF DESIRED
C
      IF(LGRAD) CALL DERIV(GEO,GRAD)
      RETURN
      END
