      SUBROUTINE COMPFG(XPARAM,INT,ESCF,FULSCF,GRAD,LGRAD)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES/NOLIST'
      DIMENSION XPARAM(MAXPAR),GRAD(MAXPAR)
      LOGICAL LGRAD, FULSCF
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR),DUMY(MAXPAR)
      COMMON /GEOSYM/ NDEP,LOCPAR(200),IDEPFN(200),LOCDEP(200)
      COMMON /GMETRY/ GEO(3,NUMATM)
      COMMON /ATHEAT/ ATHEAT
      COMMON /WMATRX/ W(N2ELEC)
      COMMON /ENUCLR/ ENUCLR
      COMMON /ELECT / ELECT
      COMMON /HMATRX/ H(MPACK)
C**************************************************************************
C
C   COMPFG CALCULATES (A) THE HEAT OF FORMATION OF THE SYSTEM, AND
C                     (B) THE GRADIENTS, IF LGRAD IS .TRUE.
C
C   ON INPUT  XPARAM = ARRAY OF PARAMETERS TO BE USED IN INTERNAL COORDS.
C             LGRAD  = .TRUE. IF GRADIENTS ARE NEEDED, .FALSE. OTHERWISE.
C             INT    = .TRUE. FOR XPARAM IN INTERNAL COORDINATES
C                      .FALSE. FOR CARTESIAN COORDINATES
C             FULSCF = .TRUE. IF FULL SCF TO BE DONE, .FALSE. OTHERWISE.
C
C   ON OUTPUT ESCF  = HEAT OF FORMATION.
C             GRAD   = ARRAY OF GRADIENTS, IF LGRAD = .TRUE.
C                      
C**************************************************************************
      COMMON /KEYWRD/KEYWRD
      CHARACTER*80 KEYWRD
      LOGICAL FIRST, DEBUG, INT, PRINT
      DIMENSION COORD(3,NUMATM)
      DATA FIRST /.TRUE./
      IF (FIRST) THEN
      FIRST=.FALSE.
      PRINT=(INDEX(KEYWRD,'COMPFG') .NE. 0)
      DEBUG=(INDEX(KEYWRD,'DEBUG') .NE. 0 .AND. PRINT)
      ENDIF
C
C SET UP COORDINATES FOR CURRENT CALCULATION
C
      IF (INT) THEN
C       PLACE THE NEW VALUES OF THE VARIABLES IN THE ARRAY GEO.
C       MAKE CHANGES IN THE GEOMETRY.
         DO 10 I=1,NVAR
            K=LOC(1,I)
            L=LOC(2,I)
10          GEO(L,K)=XPARAM(I)
C       IMPOSE THE SYMMETRY CONDITIONS + COMPUTE THE DEPENDENT-PARAMETERS,
         IF(NDEP.NE.0) CALL SYMTRY
C       NOW COMPUTE THE ATOMIC COORDINATES.
         IF( DEBUG ) THEN
             WRITE(6,FMT='('' INTERNAL COORDS'',/100(/,3F12.6))')
     +            ((GEO(J,I),J=1,3),I=1,5)
         END IF
         CALL GMETRY(GEO,COORD)
         IF( DEBUG ) THEN
             WRITE(6,FMT='('' CARTESIAN COORDS'',/100(/,3F12.6))')
     +            ((COORD(J,I),J=1,3),I=1,5)
         END IF
         CALL HCORE(COORD,H,W,ENUCLR)
      ELSE
C         LGRAD = .FALSE.
C                    DERIV DOES NOT WORK WITH COORD INPUT
      CALL HCORE(XPARAM,H,W,ENUCLR)
      END IF
C
C COMPUTE THE HEAT OF FORMATION.
C
      CALL ITER(H,W,ELECT,FULSCF)
      ESCF=(ELECT+ENUCLR)*23.061D0+ATHEAT
      IF(PRINT) WRITE(6,'(/10X,'' HEAT OF FORMATION'',F18.9)')ESCF
C
C FIND DERIVATIVES IF DESIRED
C
      IF(LGRAD) CALL DERIV(GEO,GRAD)
      RETURN 
      END
