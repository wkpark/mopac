      SUBROUTINE MAKPOL(COORD)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INCLUDE 'SIZES'
      DIMENSION COORD(3,*)
************************************************************************
*
*   MAKPOL TAKES A PRIMITIVE UNIT CELL AND GENERATES A TOTAL OF 'MERS' 
*   COPIES.  THE RESULTING GEOMETRY IS PLACED IN GEO.  ARRAYS LOC,
*   XPARAM, NA, NB, NC, SIMBOL, TXTATM, LABELS, LOCPAR, IDEPFN, AND
*   LOCDEP ARE EXPANDED TO SUIT.  ARRAY TVEC IS MODIFIED, AS ARE SCALARS
*   NVAR, NATOMS, AND NDEP.
*
*   SYMMETRY IS FORCED ON, OR ADDED ON, IN ORDER TO MAKE THE NEW MERS
*   EQUIVALENT TO THE SUPPLIED MER.
*   
************************************************************************

      CHARACTER KEYWRD*241, TXTATM*8, SIMBOL*10, LTXT*1
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOVAR/ NVAR, LOC(2,MAXPAR), IDUMY, XPARAM(MAXPAR)
      COMMON /GEOM  / GEO(3,NUMATM), XCOORD(3,NUMATM)
      COMMON /ATOMTX/ LTXT, TXTATM(NUMATM)
      COMMON /SIMBOL/ SIMBOL(MAXPAR)
      COMMON /EULER / TVEC(3,3), ID
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOSYM/ NDEP, LOCPAR(MAXPAR), IDEPFN(MAXPAR),
     1                      LOCDEP(MAXPAR)
         IOFF=0
         MERS=READA(KEYWRD,INDEX(KEYWRD,' MERS'))
         DO 270 I=1,NATOMS
  270    IF(LABELS(I).EQ.99)LABELS(I)=100
         CALL GMETRY(GEO,COORD)
         DO 280 I=1,NATOMS
  280    IF(LABELS(I).EQ.100)LABELS(I)=99
         NAN=NA(NATOMS-1)
         NBN=NB(NATOMS-1)
         NCN=NC(NATOMS-1)
         DO 330 I=2,MERS+1
            IM1=IOFF
            IOFF=IOFF+NATOMS-2
C
C   FILL THE NA, NB, AND NC ADDRESSES FOR THE NEW ATOMS
C
            DO 310 J=1,NATOMS-2
            IF(J.NE.1.AND.I.GT.MERS)GOTO 310
               SIMBOL(IOFF+J)=SIMBOL(IM1+J)
               IF(IOFF+J.NE.NATOMS-1)THEN
               NA(IOFF+J)=NA(IM1+J)+NATOMS-2
               NB(IOFF+J)=NB(IM1+J)+NATOMS-2
               NC(IOFF+J)=NC(IM1+J)+NATOMS-2
               ENDIF
               LABELS(IOFF+J)=LABELS(IM1+J)
               TXTATM(IOFF+J)=TXTATM(IM1+J)
               DO 300 K=1,3
  300          COORD(K,IOFF+J)=COORD(K,IM1+J)+TVEC(K,1)
  310       CONTINUE
            IF(I.EQ.2)THEN
C
C  SPECIAL TREATMENT FOR THE FIRST THREE ATOMS OF THE SECOND MER
C
               NA(NATOMS-1)=NAN
               NB(NATOMS-1)=NBN
               NC(NATOMS-1)=NCN
               NB(NATOMS+0)=NA(NATOMS-2)
               NC(NATOMS+0)=NB(NATOMS-2)
               NC(NATOMS+1)=NA(NATOMS-2)
            ENDIF
C#            DO 320 J=1,NATOMS-2
C#  320       WRITE(6,'(3I5,3F12.5,3I4)')I,J,LABELS(IFF+J),
C#     1(COORD(K,IOFF+J),K=1,3),
C#     2NA(IOFF+J), NB(IOFF+J), NC(IOFF+J)
  330    CONTINUE
C
C  USE ATOMS OF FIRST MER TO DEFINE THE OTHER MERS.  FOR ATOMS 1, 2, AND
C  3, USE DATA FROM THE SECOND MER.
C
         DO 360 I=1,NATOMS-2
            DO 350 K=1,3
            IF(K.GE.I)THEN
            KOFF=NATOMS-2
            JOFF=3
            ELSE
            KOFF=0
            JOFF=2
            ENDIF
               DO 340 J=JOFF,MERS+1
               IF(I.NE.1.AND.J.GT.MERS) GOTO 340
                  NDEP=NDEP+1
                  LOCPAR(NDEP)=I+KOFF
                  IDEPFN(NDEP)=K
                  LOCDEP(NDEP)=(NATOMS-2)*(J-1)+I
  340          CONTINUE
  350       CONTINUE
  360    CONTINUE
C
C   CARTESIAN COORDINATES OF THE TV
C
         LAST=(NATOMS-2)*MERS+2
         COORD(1,LAST)=COORD(1,IOFF+1)
         COORD(2,LAST)=COORD(2,IOFF+1)
         COORD(3,LAST)=COORD(3,IOFF+1)
C
C  REMOVE OPTIMIZATION FLAGS OF LAST TWO ATOMS SUPPLIED BY THE USER
C
         DO 331 I=1,6
  331    IF(LOC(1,NVAR).GT.NATOMS-2)NVAR=NVAR-1
C
C   PUT ON OPTIMIZATION FLAGES FOR FIRST THREE ATOMS OF THE SECOND MER
C
         LOC(1,NVAR+1)=NATOMS-1
         LOC(2,NVAR+1)=1
         LOC(1,NVAR+2)=NATOMS-1
         LOC(2,NVAR+2)=2
         LOC(1,NVAR+3)=NATOMS-1
         LOC(2,NVAR+3)=3
         LOC(1,NVAR+4)=NATOMS
         LOC(2,NVAR+4)=2
         LOC(1,NVAR+5)=NATOMS
         LOC(2,NVAR+5)=3
         LOC(1,NVAR+6)=NATOMS+1
         LOC(2,NVAR+6)=3
C
C  RE-DO SPECIFICATION OF THE TV
C
         LABELS(LAST-1)=99
         LABELS(LAST)=107
         TXTATM(LAST-1)=' '
         TXTATM(LAST)=' '
         NA(LAST)=1
         NB(LAST)=LAST-1
         NC(LAST)=LAST-2
         LOC(1,NVAR+7)=LAST
         LOC(2,NVAR+7)=1
C
C   CONVERT TO INTERNAL COORDINATES.  USE CONNECTIVITY CREATED HERE
C
         DEGREE=1.D0
         NA(2)=-2
         CALL XYZINT(COORD,LAST,NA,NB,NC,DEGREE,GEO)
C
C  RE-SIZE THE TRANSLATION VECTOR
C
         TVEC(1,1)=COORD(1,LAST)
         TVEC(2,1)=COORD(2,LAST)
         TVEC(3,1)=COORD(3,LAST)
C
C THE COORDINATES OF THE FIRST 3 ATOMS NEED TO BE OPTIMIZED
C
         XPARAM(NVAR+1)=GEO(1,NATOMS-1)
         XPARAM(NVAR+2)=GEO(2,NATOMS-1)
         XPARAM(NVAR+3)=GEO(3,NATOMS-1)
         XPARAM(NVAR+4)=GEO(2,NATOMS)
         XPARAM(NVAR+5)=GEO(3,NATOMS)
         XPARAM(NVAR+6)=GEO(3,NATOMS+1)
         NATOMS=LAST
         XPARAM(NVAR+7)=GEO(1,NATOMS)
         NVAR=NVAR+7
         WRITE(6,160)(I,(TVEC(J,I),J=1,3),I=1,ID)
  150    FORMAT(/,'       EXPANDED UNIT CELL TRANSLATION VECTORS',/
     1/,'              X              Y              Z')
  160    FORMAT('    T',I1,' = ',F11.7,'    ',F11.7,'    ',F11.7)
         WRITE(6,'(/,10X,A)')' EXPANDED POLYMER UNIT CELL'
         CALL GEOUT(1)
         RETURN
         END
