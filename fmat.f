      SUBROUTINE FMAT(FMATRX, NREAL, TSCF, TDER, DELDIP, HEAT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION FMATRX(*), DELDIP(3,*)
***********************************************************************
*
*  VALUE CALCULATES THE SECOND-ORDER OF THE ENERGY WITH
*        RESPECT TO THE CARTESIAN COORDINATES I AND J AND PLACES IT
*        IN FMATRX
*
*  ON INPUT NATOMS  = NUMBER OF ATOMS IN THE SYSTEM.
*           XPARAM  = INTERNAL COORDINATES OF MOLECULE STORED LINEARLY
*
*  VARIABLES USED
*           COORDL  = ARRAY OF CARTESIAN COORDINATES, STORED LINEARLY.
*           I       = INDEX OF CARTESIAN COORDINATE.
*           J       = INDEX OF CARTESIAN COORDINATE.
*
*  ON OUTPUT FMATRX = SECOND DERIVATIVE OF THE ENERGY WITH RESPECT TO
*                    CARTESIAN COORDINATES I AND J.
***********************************************************************
      COMMON /KEYWRD/ KEYWRD
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /GEOVAR/ NVAR,LOC(2,MAXPAR),IDUMY, DUMY(MAXPAR)
      COMMON /DENSTY/ P(MPACK),PDUMY(2,MPACK)
      COMMON /ATMASS/ ATMASS(NUMATM)
      COMMON /TIME  / TIME0
      COMMON /CORE  / CORE(107)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /COORD / COORD(3,NUMATM)
      COMMON /SCRACH/ EVECS(MAXPAR*MAXPAR)
      DIMENSION COLD(MAXPAR), GRAD(MAXPAR),
     1GROLD(MAXPAR), COORDL(MAXPAR), Q(NUMATM), DEL2(3), G2OLD(MAXPAR)
     2, EIGS(MAXPAR), G2RAD(MAXPAR),
     3 FCONST(MAXPAR)
      CHARACTER*80 KEYWRD
      LOGICAL DEBUG, DERIV, RESTRT, PRNT, RESFIL, ANALYT, PRECIS
      CHARACTER SPACE*1, DOTT*1, ZERO*1, NINE*1, CH*1
      EQUIVALENCE (COORD(1,1),COORDL(1))
      DATA SPACE,DOTT,ZERO,NINE /' ','.','0','9'/
      DATA FACT/6.95125D-3/
C
C    FACT IS THE CONVERSION FACTOR FROM KCAL/MOLE TO ERGS
C
C SET UP CONSTANTS AND FLAGS
      NA(1)=99
C
C  SET UP THE VARIABLES IN XPARAM ANDLOC,THESE ARE IN CARTESIAN COORDINA
C
      NUMAT=0
      DO 10 I=1,NATOMS
         IF(LABELS(I).NE.99.AND.LABELS(I).NE.107) THEN
            NUMAT=NUMAT+1
            LABELS(NUMAT)=LABELS(I)
         ENDIF
   10 CONTINUE
      NATOMS=NUMAT
C
C   THIS IS A QUICK, IF CLUMSY, WAY TO CALCULATE NUMAT, AND TO REMOVE TH
C   DUMMY ATOMS FROM THE ARRAY LABELS.
C
      NVAR=0
      DO 20 I=1,NUMAT
         DO 20 J=1,3
            NVAR=NVAR+1
            LOC(1,NVAR)=I
            LOC(2,NVAR)=J
   20 CONTINUE
      LIN=(NVAR*(NVAR+1))/2
      DO 30 I=1,LIN
   30 FMATRX(I)=0.D0
      PRNT   =(INDEX(KEYWRD,'IRC=') .EQ. 0)
      PRECIS =(INDEX(KEYWRD,'PRECIS') .NE. 0)
      ANALYT =(INDEX(KEYWRD,'ANALYT') .NE. 0)
      RESTRT =(INDEX(KEYWRD,'RESTART') .NE. 0)
      IF(INDEX(KEYWRD,'NLLSQ') .NE. 0) RESTRT=.FALSE.
      DEBUG =(INDEX(KEYWRD,'FMAT') .NE. 0)
      DERIV=(NCLOSE .EQ. NOPEN .AND. INDEX(KEYWRD,'C.I.') .EQ. 0)
      IF(PRNT)WRITE(6,'(//4X,''FIRST DERIVATIVES WILL BE USED IN THE''
     1,'' CALCULATION OF SECOND DERIVATIVES'')')
      TIME=MAXTIM
      I=INDEX(KEYWRD,' T=')
      IF(I.NE.0) THEN
         TIM=READA(KEYWRD,I)
         DO 40 J=I+3,80
            CH=KEYWRD(J:J)
            IF( CH .NE. DOTT .AND. (CH .LT. ZERO .OR. CH .GT.NINE)) THEN
               IF( CH .EQ. 'M') TIM=TIM*60
               GOTO 50
            ENDIF
   40    CONTINUE
   50    TIME=TIM
         IF(PRNT)WRITE(6,'(/10X,''TIME DEFINED FOR THIS STEP ='',F19.2,
     1    '' SECONDS'')')TIME
      ELSE
         IF(PRNT)WRITE(6,'(/10X,''DEFAULT TIME OF'',F8.2,
     1    '' SECONDS ALLOCATED FOR THIS STEP'')')TIME
      ENDIF
      TLEFT=TIME
      TLAST=TIME
      TDUMP=MAXDMP
      I=INDEX(KEYWRD,' DUMP')
      IF(I.NE.0) THEN
         TDUMP=READA(KEYWRD,I)
         DO 60 J=I+7,80
            CH=KEYWRD(J:J)
            IF( CH .NE. DOTT .AND. (CH .LT. ZERO .OR. CH .GT. NINE))
     1 THEN
               IF( CH .EQ. 'M') TDUMP=TDUMP*60
               GOTO 70
            ENDIF
   60    CONTINUE
   70    CONTINUE
      ENDIF
      RESFIL=.FALSE.
      IF(RESTRT) THEN
         DO 80 I=1,NVAR
   80    COLD(I)=COORDL(I)
         ISTART = 0
         I=0
         CALL FORSAV(TOTIME,DELDIP,ISTART,I,FMATRX, COORD, NVAR,HEAT,
     1                EVECS,JSTART,FCONST)
         KOUNTF=(ISTART*(ISTART+1))/2
         ISTART=ISTART+1
         JSTART=JSTART+1
         TIME2 = SECOND()
         IF(ISTART.GT.NVAR) GOTO 200
      ELSE
         TOTIME=0.D0
         IF (TSCF.GT.0.D0)TLEFT=TLEFT-TSCF-TDER
         ISTART=1
      ENDIF
C CALCULATE FMATRX
      IF(ISTART.GT.1) THEN
         ESTIME=(NVAR-ISTART+1)*TOTIME/(ISTART-1.D0)
      ELSE
         ESTIME=NVAR*(TSCF+TDER)*2.D0
         IF (PRECIS) ESTIME=ESTIME*2.D0
      ENDIF
      IF(TSCF.GT.0)
     1WRITE(6,'(/10X,''ESTIMATED TIME TO COMPLETE CALCULATION =''
     2,F9.2,'' SECONDS'')')ESTIME
      IF(RESTRT) THEN
         IF(ISTART.LE.NVAR)
     1    WRITE(6,'(/10X,''STARTING AGAIN AT LINE'',18X,I4)')ISTART
         WRITE(6,'(/10X,''TIME USED UP TO RESTART ='',F22.2)')TOTIME
      ENDIF
      LU=KOUNTF
      TIME1 = SECOND()
      NUMAT=NVAR/3
      DO 160 I=ISTART,NVAR
         TIME2 = SECOND()
         DELTA=1.D0/60.D0
         IF(PRECIS)THEN
C
C   DETERMINE A GOOD STEP SIZE
C
            G2OLD(1)=100.D0
            COORDL(I)=COORDL(I)+DELTA
            CALL COMPFG(COORDL, .TRUE., ESCF, .TRUE., G2OLD, .TRUE.)
            COORDL(I)=COORDL(I)-DELTA
            DELTA=DELTA*10.D0/SQRT(DOT(G2OLD,G2OLD,NVAR))
C#         WRITE(6,'(A,F12.5)')' DELTA:',DELTA
            G2OLD(1)=100.D0
            COORDL(I)=COORDL(I)+DELTA
            CALL COMPFG(COORDL, .TRUE., ESCF, .TRUE., G2OLD, .TRUE.)
C#         WRITE(6,*)' GNORM:',SQRT(DOT(G2OLD,G2OLD,NVAR))
            COORDL(I)=COORDL(I)-DELTA*2.D0
            G2RAD(1)=100.D0
            CALL COMPFG(COORDL, .TRUE., HEATAA, .TRUE., G2RAD, .TRUE.)
            COORDL(I)=COORDL(I)+DELTA
         ENDIF
         COORDL(I)=COORDL(I)+0.5D0*DELTA
         GROLD(1)=100.D0
         CALL COMPFG(COORDL, .TRUE., ESCF, .TRUE., GROLD, .TRUE.)
C#         WRITE(6,*)' GNORM:',SQRT(DOT(GROLD,GROLD,NVAR))
         CALL CHRGE(P,Q)
         DO 90 II=1,NUMAT
   90    Q(II)=CORE(LABELS(II))-Q(II)
         SUM = DIPOLE(P,Q,COORDL,DELDIP(1,I))
         COORDL(I)=COORDL(I)-DELTA
         GRAD(1)=100.D0
         CALL COMPFG(COORDL, .TRUE., HEATAA, .TRUE., GRAD, .TRUE.)
         COORDL(I)=COORDL(I)+DELTA*0.5D0
         CALL CHRGE(P,Q)
         DO 100 II=1,NUMAT
  100    Q(II)=CORE(LABELS(II))-Q(II)
         SUM = DIPOLE(P,Q,COORDL,DEL2)
         DO 110 II=1,3
  110    DELDIP(II,I)=(DELDIP(II,I)-DEL2(II))*0.5D0/DELTA
         LL=LU+1
         LU=LL+I-1
         L=0
         IF(PRECIS)THEN
            DO 120 KOUNTF=LL,LU
               L=L+1
               FMATRX(KOUNTF)=FMATRX(KOUNTF)+
     1         (8.D0*(GROLD(L)-GRAD(L))-(G2OLD(L)-G2RAD(L)))
     2          *0.25D0/DELTA*FACT/6.D0
  120       CONTINUE
            L=L-1
            DO 130 K=I,NVAR
               L=L+1
               KK=(K*(K-1))/2+I
               FMATRX(KK)=FMATRX(KK)+
     1         (8.D0*(GROLD(L)-GRAD(L))-(G2OLD(L)-G2RAD(L)))
     2          *0.25D0/DELTA*FACT/6.D0
  130       CONTINUE
         ELSE
            DO 140 KOUNTF=LL,LU
               L=L+1
               FMATRX(KOUNTF)=FMATRX(KOUNTF)+
     1         ((GROLD(L)-GRAD(L)))
     2          *0.25D0/DELTA*FACT
  140       CONTINUE
            L=L-1
            DO 150 K=I,NVAR
               L=L+1
               KK=(K*(K-1))/2+I
               FMATRX(KK)=FMATRX(KK)+
     1         ((GROLD(L)-GRAD(L)))
     2          *0.25D0/DELTA*FACT
  150       CONTINUE
         ENDIF
         TIME3 = SECOND()
         TSTEP=TIME3-TIME2
         TOTIME= TOTIME+TSTEP
         TLEFT= TLEFT-TSTEP
         IF(RESFIL)THEN
            WRITE(6,'('' STEP:'',I4,'' RESTART FILE WRITTEN, INTEGRAL ='
     1',F10.2,'' TIME LEFT:'',F10.2)')I,TOTIME,TLEFT
            RESFIL=.FALSE.
         ELSE
            WRITE(6,'('' STEP:'',I4,'' TIME ='',F9.2,'' SECS, INTEGRAL =
     1'',F10.2,'' TIME LEFT:'',F10.2)')I,TSTEP,TOTIME,TLEFT
         ENDIF
         IF(DERIV) THEN
            ESTIM = TOTIME/I
         ELSE
            ESTIM = TOTIME*2.D0/I
         ENDIF
         IF(TLAST-TLEFT.GT.TDUMP)THEN
            TLAST=TLEFT
            RESFIL=.TRUE.
            JSTART=1
            II=I
            CALL FORSAV(TOTIME,DELDIP,II,NVAR,FMATRX, COORD,NVAR,HEAT,
     1                EVECS,JSTART,FCONST)
         ENDIF
         IF(I.NE.NVAR.AND.TLEFT-10.D0 .LT. ESTIM) THEN
            WRITE(6,'(//10X,''- - - - - - - TIME UP - - - - - - -'',//)'
     1)
            WRITE(6,'(/10X,'' POINT REACHED ='',I4)')I
            WRITE(6,'(/10X,'' RESTART USING KEY-WORD "RESTART"'')')
            WRITE(6,'(10X,''ESTIMATED TIME FOR THE NEXT STEP ='',F8.2,
     1'' SECONDS'')')ESTIM
            JSTART=1
            II=I
            CALL FORSAV(TOTIME,DELDIP,II,NVAR,FMATRX, COORD,NVAR,HEAT,
     1                EVECS,JSTART,FCONST)
            WRITE(6,'(//10X,''FORCE MATRIX WRITTEN TO DISK'')')
            STOP
         ENDIF
  160 CONTINUE
C#      CALL FORSAV(TOTIME,DELDIP,NVAR,NVAR,FMATRX, COORD,NVAR,HEAT,
C#     +                EVECS,JSTART,FCONST)
      IF(DERIV) GOTO 290
      WRITE(6,'(//10X,'' STARTING TO CALCULATE FORCE CONSTANTS'',/)')
      CALL FRAME(FMATRX,NUMAT,0,EIGS)
      CALL RSP(FMATRX,NVAR,NVAR,EIGS,EVECS)
      IF(DEBUG) THEN
         WRITE(6,'(''   EIGENVECTORS FROM FIRST CALCULATION'')')
         CALL MATOUT(EVECS,EIGS,NVAR,NVAR,NVAR)
      ENDIF
      L=0
      DO 180 I=1,NVAR
         DO 180 J=1,I
            L=L+1
            SUM=0.D0
            DO 170 K=1,NREAL
               K1=(K-1)*NVAR+I
               K2=(K-1)*NVAR+J
  170       SUM=SUM+EVECS(K1)*EIGS(K)*EVECS(K2)
  180 FMATRX(L)=SUM
      CALL FRAME(FMATRX,NUMAT,0,EIGS)
      CALL RSP(FMATRX,NVAR,NVAR,EIGS,EVECS)
C#      CALL MATOUT(EVECS,EIGS,NVAR,NVAR,NVAR)
      JSTART=1
      DO 190 I=1,NVAR
  190 COLD(I)=COORDL(I)
  200 IF(DERIV) GOTO 290
      DELTA=0.025D0
      L=(JSTART-1)*NVAR
      DO 250 ILOOP=JSTART,NVAR
         J=L
         DO 210 I=1,NVAR
            J=J+1
  210    COORDL(I)=COLD(I)+EVECS(J)*DELTA
         CALL COMPFG(COORDL, .TRUE., HEATA, .TRUE., GRAD, .FALSE.)
         HEATA=HEATA-HEAT
         J=L
         DO 220 I=1,NVAR
            J=J+1
  220    COORDL(I)=COLD(I)-EVECS(J)*DELTA
         CALL COMPFG(COORDL, .TRUE., HEATB, .TRUE., GRAD, .FALSE.)
         HEATB=HEATB-HEAT
         J=L
         DO 230 I=1,NVAR
            J=J+1
  230    COORDL(I)=COLD(I)+EVECS(J)*DELTA*2
         CALL COMPFG(COORDL, .TRUE., HEATAA, .TRUE., GRAD, .FALSE.)
         HEATAA=HEATAA-HEAT
         J=L
         DO 240 I=1,NVAR
            J=J+1
  240    COORDL(I)=COLD(I)-EVECS(J)*DELTA*2
         CALL COMPFG(COORDL, .TRUE., HEATBB, .TRUE., GRAD, .FALSE.)
         HEATBB=HEATBB-HEAT
         SUM=( (HEATA+HEATB)*16 - (HEATAA+HEATBB) )/12.D0
     1/DELTA*FACT/DELTA*0.5D0
         FCONST(ILOOP)=SUM*0.5D0
         L=L+NVAR
         TIME3 = SECOND()
         TSTEP=TIME3-TIME2
         TIME2=TIME3
         TOTIME= TOTIME+TSTEP
         TLEFT= TLEFT-TSTEP
         WRITE(6,'('' STEP:'',I4,'' TIME ='',F9.2,'' SECS, INTEGRAL ='',
     1F10.2,'' TIME LEFT:'',F10.2)')ILOOP,TSTEP,TOTIME,TLEFT
         ESTIM = TSTEP*5.D0
C
C    5.0 IS A SAFETY FACTOR
C
         IF(TLAST-TLEFT.GT.TDUMP)THEN
            TLAST=TLEFT
            RESFIL=.TRUE.
            IFOR=ILOOP
            IX=NVAR+2
*
* VALUE OF IX IS NOT IMPORTANT. SHOULD NOT BE 0 OR NVAR
*
            CALL FORSAV(TOTIME,DELDIP,IX,NVAR,FMATRX, COORD,NVAR,HEAT,
     1                EVECS,IFOR,FCONST)
         ENDIF
         IF(ILOOP.NE.NVAR.AND.TLEFT-10.D0 .LT. ESTIM) THEN
            WRITE(6,'(//10X,''- - - - -  TIME  LIMIT - - - - -'')')
            WRITE(6,'(/10X,'' POINT REACHED ='',I4)')ILOOP
            WRITE(6,'(/10X,'' RESTART USING KEY-WORD "RESTART"'')')
            WRITE(6,'(10X,''ESTIMATED TIME FOR THE NEXT STEP ='',F8.2,
     1'' SECONDS'')')ESTIM
            IFOR=ILOOP
            IX=NVAR+2
*
* VALUE OF IX IS NOT IMPORTANT. SHOULD NOT BE 0 OR NVAR
*
            CALL FORSAV(TOTIME,DELDIP,IX,NVAR,FMATRX, COORD,NVAR,HEAT,
     1                EVECS,IFOR,FCONST)
         ENDIF
  250 CONTINUE
      L=0
      DO 270 I=1,NVAR
         DO 270 J=1,I
            L=L+1
            SUM=0.D0
            DO 260 K=1,NVAR
               K1=(K-1)*NVAR+I
               K2=(K-1)*NVAR+J
  260       SUM=SUM+EVECS(K1)*FCONST(K)*EVECS(K2)
  270 FMATRX(L)=SUM*2.D0
      DO 280 I=1,NVAR
  280 COORDL(I)=COLD(I)
  290 CONTINUE
      DO 300 I=1,NUMAT
      IF(ATMASS(I).LT.1.D-20)THEN
      CALL FORSAV(TOTIME,DELDIP,NVAR,NVAR,FMATRX, COORD,NVAR,HEAT,
     2                EVECS,ILOOP,FCONST)
      WRITE(6,'(A)')' AT LEAST ONE ATOM HAS A ZERO MASS. A RESTART'
      WRITE(6,'(A)')' FILE HAS BEEN WRITTEN AND THE JOB STOPPED'
      STOP
      ENDIF
  300 CONTINUE
      IF(ISTART.LE.NVAR .AND. INDEX(KEYWRD,'ISOTOPE') .NE. 0)
     1CALL FORSAV(TOTIME,DELDIP,NVAR,NVAR,FMATRX, COORD,NVAR,HEAT,
     2                EVECS,ILOOP,FCONST)
      RETURN
      END
