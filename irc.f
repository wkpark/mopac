      SUBROUTINE IRC
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
************************************************************************
*                                                                      *
*    IRC IS DESIGNED TO FOLLOW A REACTION PATH FROM THE TRANSITION     *
*    STATE.  TWO MODES ARE SUPPORTED, FIRST: GAS PHASE:- AS THE SYSTEM *
*    MOVES FROM THE T/S THE MOMENTUM OF THE ATOMS IS STORED AND THE    *
*    POSITION OF THE ATOMS IS RELATED TO THE OLD POSITION BY (A) THE   *
*    CURRENT VELOCITY OF THE ATOM, AND (B) THE FORCES ACTING ON THAT   *
*    ATOM.  THE SECOND MODE IS CONDENSED PHASE, IN WHICH THE ATOMS MOVE*
*    IN RESPONSE TO THE FORCES ACTING ON THEM. I.E. INFINITELY DAMPED  *
*                                                                      *
************************************************************************
      INCLUDE 'SIZES'
      COMMON /KEYWRD/ KEYWRD
      COMMON /TITLES/ KOMENT,TITLE
      COMMON /GRADNT/ GRAD(MAXPAR),GNORM
      COMMON /GEOSYM/ NDEP,LOCPAR(MAXPAR),IDEPFN(MAXPAR),LOCDEP(MAXPAR)
      COMMON /ELEMTS/ ELEMNT(107)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /ATMASS/ ATMASS(NUMATM)
      COMMON /GEOVAR/ NVAR, LOC(2,MAXPAR), IDUMY, XPARAM(MAXPAR)
      COMMON /GEOKST/ NATOMS,LABELS(NUMATM),
     1                NA(NUMATM),NB(NUMATM),NC(NUMATM)
      COMMON /CORE  / CORE(107)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /DENSTY/ P(MPACK),PA(MPACK),PB(MPACK)
      DIMENSION Q(3), COORDS(3,NUMATM), IEL1(3), Q2(NUMATM)
      CHARACTER KEYWRD*80, KOMENT*80, TITLE*80, ELEMNT*2, ALPHA*2
      CHARACTER SPACE*1, CHDOT*1, ZERO*1, NINE*1, CH*1
      DIMENSION VELO0(MAXPAR), VELO1(MAXPAR), VELO2(MAXPAR),
     1VELO3(MAXPAR),
     2VOLD(MAXPAR), COORD(3,NUMATM), GROLD2(MAXPAR),
     3GROLD(MAXPAR), PAROLD(MAXPAR), PARREF(MAXPAR),
     4 SQRTMS(MAXPAR)
      LOGICAL INT, ADDK
      EQUIVALENCE (COORD,VOLD)
      DATA VELO0/MAXPAR*0.D0/, INT/.TRUE./, VOLD/MAXPAR*0.D0/
      DATA ADDK/.TRUE./
      DATA SPACE,CHDOT,ZERO,NINE /' ','.','0','9'/
      CLOSE (5)
      TNOW=SECOND()
      IF(INDEX(KEYWRD,' PREC').NE.0)THEN
         ACCU=1.D0
      ELSE
         ACCU=0.1D0
      ENDIF
      LPOINT=0
      STEPT=0.D0
      STEPH=0.D0
      STEPX=0.D0
      IF(INDEX(KEYWRD,' T-PRIO').NE.0)THEN
         IF(INDEX(KEYWRD,' T-PRIORITY=').NE.0)THEN
            STEPT=READA(KEYWRD,INDEX(KEYWRD,'T-PRIO')+5)
         ELSE
            STEPT=0.1D0
         ENDIF
         WRITE(6,'(/,'' TIME PRIORITY, INTERVAL ='',F4.1,
     1'' FEMTOSECONDS'',/)')STEPT
         STEPT=STEPT*1.D-15
      ELSEIF(INDEX(KEYWRD,' H-PRIO').NE.0)THEN
         IF(INDEX(KEYWRD,' H-PRIORITY=').NE.0)THEN
            STEPH=READA(KEYWRD,INDEX(KEYWRD,'H-PRIO')+5)
         ELSE
            STEPH=0.1D0
         ENDIF
         WRITE(6,'(/,'' KINETIC ENERGY PRIORITY, STEP ='',F4.1,
     1'' KCAL/MOLE'',/)')STEPH
      ELSEIF(INDEX(KEYWRD,' X-PRIO').NE.0)THEN
         IF(INDEX(KEYWRD,' X-PRIORITY=').NE.0)THEN
            STEPX=READA(KEYWRD,INDEX(KEYWRD,'X-PRIO')+5)
         ELSE
            STEPX=0.05D0
         ENDIF
         WRITE(6,'(/,'' GEOMETRY PRIORITY, STEP ='',F5.2,
     1'' ANGSTROMS'',/)')STEPX
      ENDIF
      IF(INDEX(KEYWRD,' SYMME').NE.0)THEN
         WRITE(6,*)'  SYMMETRY SPECIFIED, BUT CANNOT BE USED IN IRC'
         NDEP=0
      ENDIF
      IF(INDEX(KEYWRD,' XYZ').EQ.0)THEN
         CALL GMETRY(GEO,COORD)
         L=0
         DO 10 I=1,NUMAT
            LABELS(I)=NAT(I)
            SUM=SQRT(ATMASS(NAT(I)))
            DO 10 J=1,3
               L=L+1
               SQRTMS(L)=SUM
   10    GEO(J,I)=COORD(J,I)
         NA(1)=99
      ENDIF
      L=0
      DO 20 I=1,NUMAT
         DO 20 J=1,3
            L=L+1
            LOC(1,L)=I
            LOC(2,L)=J
   20 XPARAM(L)=GEO(J,I)
      NVAR=NUMAT*3
C
C DETERMINE DAMPING FACTOR
C
      IF(INDEX(KEYWRD,'IRC=').NE.0) THEN
         HALF=READA(KEYWRD,INDEX(KEYWRD,'IRC='))
         WRITE(6,'(//10X,'' DAMPING FACTOR FOR KINETIC ENERGY ='',F12.6)
     1')HALF
      ELSE
         HALF=1.D6
      ENDIF
      HALF=MAX(0.0001D0,HALF)
C
C DETERMINE EXCESS KINETIC ENERGY
C
      IF(INDEX(KEYWRD,'KINE').NE.0) THEN
         ADDONK=READA(KEYWRD,INDEX(KEYWRD,'KINE'))
         WRITE(6,'(//10X,'' EXCESS KINETIC ENERGY ENTERED INTO SYSTEM ='
     1',F12.6)')ADDONK
      ELSE
         ADDONK=0.D0
      ENDIF
C
C   LOOP OVER TIME-INTERVALS OF DELTAT SECOND
C
      DELTAT=1.D-16
      QUADR=1.D0
      TOTIME=0.D0
      ETOT=0.D0
      ESCF=0.D0
      CONST=1.D0
         I=INDEX(KEYWRD,' T=')
         IF(I.NE.0) THEN
            TIM=READA(KEYWRD,I)
            DO 11 J=I+3,80
               CH=KEYWRD(J:J)
               IF( CH .NE. CHDOT .AND. (CH .LT. ZERO .OR. CH .GT. NINE))
     1 THEN
                  IF( CH .EQ. 'M') TIM=TIM*60
                  GOTO 21
               ENDIF
   11       CONTINUE
C             4 SECONDS TO LOAD IN EXECUTABLE!
   21       TLEFT=TIM-4
         ELSE
         TLEFT=3596
         ENDIF
      IF( INDEX(KEYWRD,'REST').NE.0)THEN
C
C  RESTART FROM A PREVIOUS RUN
C
      OPEN(UNIT=9,FILE='FOR009',STATUS='UNKNOWN',FORM='FORMATTED')
      REWIND 9
      READ(9,'(A)')ALPHA
      READ(9,'(8F10.6)')(XPARAM(I),I=1,NVAR)
      READ(9,'(A)')ALPHA
      READ(9,'(8F10.1)')(VELO0(I),I=1,NVAR)
      READ(9,'(A)')ALPHA
      READ(9,*)(GRAD(I),I=1,NVAR)
      READ(9,*)(GROLD(I),I=1,NVAR)
      READ(9,*)(GROLD2(I),I=1,NVAR)
      READ(9,*)(PARREF(I),I=1,NVAR)
      READ(9,*)ETOT,ESCF,EKIN,DELOLD,DELTAT,DLOLD2,ILOOP,
     +TOTIME,OLDT,TREF,OLDH,OLDT,REFSCF
      WRITE(6,'(//10X,''CALCULATION RESTARTED, CURRENT'',
     +'' KINETIC ENERGY='',F10.5,//)')EKIN
      ELSE
      ILOOP=1
      ENDIF
      IUPPER=ILOOP+399
      DO 230 ILOOP=ILOOP,IUPPER
C
C  MOVEMENT OF ATOMS WILL BE PROPORTIONAL TO THE AVERAGE VELOCITIES
C  OF THE ATOMS BEFORE AND AFTER TIME INTERVAL
C
         DO 30 I=1,NVAR
   30    VOLD(I)=VELO0(I)
C
C   KINETIC ENERGY = 1/2 * M * V * V
C                  = 0.5 / (4.184D10) * M * V * V
C   NEW VELOCITY = OLD VELOCITY + GRADIENT * TIME / MASS
C                = KCAL/ANGSTROM*SECOND/(ATOMIC WEIGHT)
C                =4.184*10**10(ERGS)*10**8(PER CM)*DELTAT(SECONDS)
C   NEW POSITION = OLD POSITION - AVERAGE VELOCITY * TIME INTERVAL
C
C
C   ESTABLISH REFERENCE TOTAL ENERGY
C
         ERROR=(ETOT-(EKIN+ESCF))
         IF(ILOOP.GT.2)THEN
            QUADR = 1.D0+ERROR/(EKIN*CONST+0.001D0)*0.5D0
C#      WRITE(6,'(''   QUADR'',3F13.6)')QUADR,ERROR,EKIN
            QUADR = MIN(1.3D0,MAX(0.8D0,QUADR))
         ELSE
            QUADR=1.D0
         ENDIF
         IF(EKIN.GT.0.2.AND.ADDK)THEN
            ETOT=ETOT+ADDONK
            ADDK=.FALSE.
         ENDIF
         EKOLD=EKIN
         IF(EKIN.GT.0.2.AND.ADDK)THEN
C
C   DUMP IN EXCESS KINETIC ENERGY
C
            ETOT=ETOT+ADDONK
            ADDK=.FALSE.
         ENDIF
         DLOLD2=DELOLD
         DELOLD=DELTAT
C
C  CALCULATE THE DURATION OF THE NEXT STEP.
C
         IF(CONST.GT.0.99D0) THEN
            DELTAT= MAX(DELTAT*(MIN(1.5D0,0.0005D0*ACCU/ABS(ERROR+1.D-11
     1))),1.D-16)
         ELSE
            DELTAT=ACCU*1.D-15
         ENDIF
         EROLD=ERROR
C
C  IF DAMPING IS USED, CALCULATE THE NEW TOTAL ENERGY AND
C  THE RATIO FOR REDUCING THE KINETIC ENERGY
C
         CONST=0.5D0**(DELTAT*1.D14/HALF)
         ETOT=ETOT-EKIN*(1-CONST)
         CONST=SQRT(CONST)
C
         VELVEC=0.D0
         EKIN=0.D0
C#      WRITE(6,*)' VELOCITY, FIRST DERIV AND SECOND DERIV'
C#      WRITE(6,'(6E13.4)')(VELO0(J),J=1,NVAR)
         DELTA1=DELTAT+DELOLD
         DELTA2=DELTA1+DLOLD2
         DO 40 I=1,NVAR
C
C   CALCULATE COMPONENTS OF VELOCITY AS
C   V = V(0) + V'*T + V"*T*T
C   WE NEED ALL THREE TERMS, V(0), V' AND V"
C
            VELO1(I) = 1.D0/ATMASS(LOC(1,I))*GRAD(I)
            VELO2(I) = 1.D0/ATMASS(LOC(1,I))*
     1                 (GRAD(I)-GROLD(I))/DELOLD
            IF(ILOOP.GT.3) VELO3(I) = 2.D0/ATMASS(LOC(1,I))*
     1      ((DELTAT-DELTA2)*(GRAD(I)-GROLD(I))
     2      -(DELTAT-DELTA1)*(GRAD(I)-GROLD2(I)))/
     3      ((1.D30*DELTAT**2-1.D30*DELTA1**2)*(DELTAT-DELTA2)
     4      -(1.D30*DELTAT**2-1.D30*DELTA2**2)*(DELTAT-DELTA1))
C
C  MOVE ATOMS THROUGH DISTANCE EQUAL TO VELOCITY * DELTA-TIME, NOTE
C  VELOCITY CHANGES FROM START TO FINISH, THEREFORE AVERAGE.
C
            PAROLD(I)=XPARAM(I)
            XPARAM(I)=XPARAM(I)
     1             -1.D8*(DELTAT*VELO0(I)
     2             +0.5D0*DELTAT**2*VELO1(I)
     3             +0.3333333D0*DELTAT**3*VELO2(I)
     4             +0.25D0*DELTAT**2*(1.D30*DELTAT**2)*VELO3(I))
C
C   CORRECT ERRORS DUE TO CUBIC COMPONENTS IN ENERGY GRADIENT,
C   ALSO TO ADD ON EXCESS ENERGY, IF NECESSARY.
C
            VELO0(I)=VELO0(I) *  QUADR
            VELVEC=VELVEC+VELO0(I)**2
C
C   MODIFY VELOCITY IN LIGHT OF CURRENT ENERGY GRADIENTS.
C
C   VELOCITY = OLD VELOCITY + (DELTA-T / ATOMIC MASS) * CURRENT GRADIENT
C                           + 1/2 *(DELTA-T * DELTA-T /ATOMIC MASS) *
C                             (SLOPE OF GRADIENT)
C              SLOPE OF GRADIENT = (GRAD(I)-GROLD(I))/DELOLD
C
C
C   THIS EXPRESSION IS ACCURATE TO SECOND ORDER IN TIME.
C
            VELO0(I) = VELO0(I) + DELTAT*VELO1(I) + 0.5D0*DELTAT**2*VELO
     12(I)                 + 0.3333333D0*DELTAT*(1.D30*DELTAT**2)*VELO3(
     2I)
C
C  CALCULATE KINETIC ENERGY (IN 2*ERGS AT THIS POINT)
C
            EKIN=EKIN+VELO0(I)**2*ATMASS(LOC(1,I))
C#      WRITE(6,*)'   VELVEC:',VELVEC,' QUADR:',QUADR
C#      WRITE(6,*)'  VELOCITY DUE TO GRADIENT TERM'
C#      WRITE(6,'(6F13.1)')(VELO0(I)-VOLD(I),I=1,NVAR)
C#      WRITE(6,*)'  VELOCITY AFTER ADDING IN GRADIENT TERM'
C#      WRITE(6,'(6F13.1)')(VELO0(I),I=1,NVAR)
   40    CONTINUE
C#      WRITE(6,'(6E13.4)')(VELO1(J),J=1,NVAR)
C#      WRITE(6,'(6E13.4)')(VELO2(J),J=1,NVAR)
C#      WRITE(6,'(6E13.4)')(VELO3(J),J=1,NVAR)
C
C  CONVERT ENERGY INTO KCAL/MOLE
C
         EKIN=0.5*EKIN/4.184D10
C
C STORE OLD GRADIENTS FOR DELTA - VELOCITY CALCULATION
C
         DO 50 I=1,NVAR
            GROLD2(I)=GROLD(I)
            GROLD(I)=GRAD(I)
   50    GRAD(I)=0.D0
C
C   CALCULATE ENERGY AND GRADIENTS
C
         CALL COMPFG(XPARAM,INT,ESCF,.TRUE.,GRAD,.TRUE.)
C
C   CONVERT GRADIENTS INTO ERGS/CM
C
         DO 60 I=1,NVAR
   60    GRAD(I)=GRAD(I)*4.184D18
C
C   SPECIAL TREATMENT FOR FIRST POINT - SET "OLD" GRADIENTS EQUAL TO
C   CURRENT GRADIENTS.
C
         IF(ILOOP.EQ.1) THEN
            DO 70 I=1,NVAR
   70       GROLD(I)=GRAD(I)
         ENDIF
C
C   GO THROUGH THE CRITERIA FOR DECIDING WHETHER OR NOT TO PRINT THIS
C   POINT.  IF YES, THEN ALSO CALCULATE THE EXACT POINT AS A FRACTION
C   BETWEEN THE LAST POINT AND THE CURRENT POINT
C
         FRACT=-1.D-4
         NFRACT=1
         IF(OLDH.EQ.0)OLDH=ESCF
         IF(STEPH.EQ.0)GOTO 80
C
C   CRITERION FOR PRINTING RESULTS  IS A CHANGE IN HEAT OF FORMATION =
C   -CHANGE IN KINETIC ENERGY
C
         IF(REFSCF.EQ.0.D0) THEN
            I=ESCF/STEPH
            REFSCF=I*STEPH
         ENDIF
         IF(ABS(ESCF-REFSCF).GT.STEPH)THEN
C
C   FRACT IS THE FRACTIONAL DISTANCE FROM THE OLD TO THE NEW POINT
C   FOR THE FIRST POINT TO BE PRINTED
C
C   FINCR IS THE DISTANCE BETWEEN ANY TWO POINTS TO BE PRINTED
C
C   NFRACT IS THE NUMBER OF POINTS TO BE PRINTED IN THE CURRENT DOMAIN
C
            FRACT=(REFSCF+SIGN(STEPH,ESCF-OLDH)-OLDH)/(ESCF-OLDH)
            FINCR=STEPH/ABS(ESCF-OLDH)
            NFRACT=(ESCF-REFSCF)/STEPH
            REFSCF=REFSCF+STEPH*NFRACT
            NFRACT=ABS(NFRACT)
         ENDIF
         GOTO 140
   80    IF(STEPT.EQ.0.D0) GOTO 90
C
C   CRITERION FOR PRINTING RESULTS IS A CHANGE IN TIME.
C
         IF(ABS(TOTIME-TREF).GT.STEPT)THEN
            FRACT = (TREF+STEPT-OLDT) / (TOTIME-OLDT)
            FINCR=STEPT/(TOTIME-OLDT)
            NFRACT= (TOTIME-TREF)/STEPT
            TREF=TREF+NFRACT*STEPT
         ENDIF
         GOTO 140
   90    IF(STEPX.EQ.0.D0)GOTO 130
C
C   CRITERION FOR PRINTING RESULTS IS A CHANGE IN GEOMETRY.
C
         XOLD=XNOW
         XNOW=0.D0
         L=0
         DO 110 I=1,NUMAT
            SUM=0.D0
            DO 100 J=1,3
               L=L+1
  100       SUM=SUM+(PARREF(L)-XPARAM(L))**2
  110    XNOW=XNOW+SQRT(SUM)
         IF(XNOW.GT.STEPX)THEN
            NFRACT=XNOW/STEPX
            SHIFT = STEPX*NFRACT
            DO 120 J=1,NVAR
  120       PARREF(J)=PARREF(J)+(XPARAM(J)-PARREF(J))/XNOW*NFRACT*STEPX
            FRACT = (STEPX-XOLD)/(XNOW-XOLD)
            FINCR=STEPX/(XNOW-XOLD)
            XNOW=XNOW-NFRACT*STEPX
            IF(ILOOP.EQ.1)NFRACT=1
         ENDIF
         GOTO 140
C
C   PRINT EVERY POINT.
C
  130    FRACT=1.0
  140    IF(ILOOP.NE.1.AND.FRACT.LT.0.D0)GOTO 210
         FRACT=FRACT-FINCR
C
C  LOOP OVER ALL POINTS IN CURRENT DOMAIN
C
         DO 220 II=1,NFRACT
            FRACT=FRACT+FINCR
            TPOINT = (1.D0-FRACT)*OLDT + FRACT*TOTIME
            HPOINT = (1.D0-FRACT)*OLDH + FRACT*ESCF
            POINTK = (1.D0-FRACT)*OLDK + FRACT*EKIN
            L=0
            DO 150 I=1,NUMAT
               DO 150 J=1,3
                  L=L+1
  150       GEO(J,I)=(1.D0-FRACT)*PAROLD(L) + FRACT*XPARAM(L)
            LPOINT=LPOINT+1
            IF(ETOT.EQ.0)ETOT=ESCF+EKIN
            IF(LPOINT.EQ.1)THEN
               WRITE(6,'(//,'' TIME IN FEMTOSECONDS  POINT  POTENTIAL +'
     1','' KINETIC  =  TOTAL     ERROR    REF%'')')
            ELSE
               WRITE(6,'(//,'' TIME IN FEMTOSECONDS  POINT  POTENTIAL +'
     1','' KINETIC  =  TOTAL     ERROR    REF'')')
            ENDIF
            WRITE(6,'(F10.3,I16,F12.5,F11.5,F11.5, F10.5,''   ''
     1,I3,''%'')')TPOINT*1.D15, ILOOP, HPOINT, POINTK, HPOINT+POINTK,
     2HPOINT+POINTK-ETOT, LPOINT
C#     1,HPOINT+POINTK-ETOT,HPOINT,POINTK,HPOINT+POINTK
            WRITE(6,*)'                CARTESIAN GEOMETRY
     1 '//'VELOCITY (IN CM/SEC)'
            WRITE(6,*)'  ATOM        X          Y          Z
     1    '//'X          Y          Z'
            DO 160 I=1,NUMAT
               LL=(I-1)*3+1
               LU=LL+2
               HPOINT = (1.D0-FRACT)*OLDH + FRACT*ESCF
               WRITE(6,'(I4,3X,A2,3F11.5,2X,3F11.1)')
     1I, ELEMNT(NAT(I)),(GEO(J,I),J=1,3),
     2((1.D0-FRACT)*VOLD(L) + FRACT*VELO0(L),L=LL,LU)
  160       CONTINUE
            CALL CHRGE(P,Q2)
            DO 170 I=1,NUMAT
               L=NAT(I)
  170       Q2(I)=CORE(L) - Q2(I)
            CALL XYZINT(GEO,NATOMS,NA,NB,NC,1.D0,COORDS)
            DEGREE=57.29577951D0
            COORDS(2,1)=0.D0
            COORDS(3,1)=0.D0
            COORDS(1,1)=0.D0
            COORDS(2,2)=0.D0
            COORDS(3,2)=0.D0
            COORDS(3,3)=0.D0
            IVAR=1
            NA(1)=0
            L=0
            WRITE(6,'(//10X,''FINAL GEOMETRY OBTAINED'',33X,''CHARGE'')'
     1)
            WRITE(6,'(1X,A)')KEYWRD,KOMENT,TITLE
            DO 200 I=1,NATOMS
               J=I/26
               ALPHA(1:1)=CHAR(ICHAR('A')+J)
               J=I-J*26
               ALPHA(2:2)=CHAR(ICHAR('A')+J-1)
               DO 180 J=1,3
  180          IEL1(J)=0
  190          CONTINUE
               IF(LOC(1,IVAR).EQ.I) THEN
                  IEL1(LOC(2,IVAR))=1
                  IVAR=IVAR+1
                  GOTO 190
               ENDIF
               IF(I.LT.4) THEN
                  IEL1(3)=0
                  IF(I.LT.3) THEN
                     IEL1(2)=0
                     IF(I.LT.2) THEN
                        IEL1(1)=0
                     ENDIF
                  ENDIF
               ENDIF
               IF(I.EQ.LATOM)IEL1(LPARAM)=-1
               Q(1)=COORDS(1,I)
               Q(2)=COORDS(2,I)*DEGREE
               Q(3)=COORDS(3,I)*DEGREE
               IF(LABELS(I).LT.99)THEN
                  L=L+1
                  WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3,F13.4,I5,A)')
     1    ELEMNT(LABELS(I)),(Q(K),IEL1(K),K=1,3),NA(I),NB(I),NC(I),Q2(L)
     2    ,LPOINT,ALPHA//'*'
               ELSE
                  WRITE(6,'(2X,A2,3(F12.6,I3),I4,2I3,13X,I5,A)')
     1    ELEMNT(LABELS(I)),(Q(K),IEL1(K),K=1,3),NA(I),NB(I),NC(I)
     2    ,LPOINT,ALPHA//'%'
               ENDIF
  200       CONTINUE
            NA(1)=99
            WRITE(6,*)
            WRITE(6,*)
  210       CONTINUE
  220    CONTINUE
         OLDT=TOTIME
         TOTIME=TOTIME+DELTAT
         OLDH=ESCF
         XOLD=XNOW
         OLDK=EKIN
      OLDTIM=TNOW
      TNOW=SECOND()
      TCYCLE=TNOW-OLDTIM
      TLEFT=TLEFT-TCYCLE      
      IF (ILOOP.EQ.IUPPER.OR.TLEFT.LT.3*TCYCLE) THEN
      OPEN(UNIT=9,FILE='FOR009',STATUS='UNKNOWN',FORM='FORMATTED')
      REWIND 9
      OPEN(UNIT=10,FILE='FOR010',STATUS='UNKNOWN',FORM='UNFORMATTED')
      REWIND 10
      WRITE(9,'(A)')' CARTESIAN GEOMETRY PARAMETERS IN ANGSTROMS'
      WRITE(9,'(8F10.6)')(XPARAM(I),I=1,NVAR)
      WRITE(9,'(A)')' VELOCITY FOR EACH CARTESIAN COORDINATE, IN CM/SEC'
      WRITE(9,'(8F10.1)')(VELO0(I),I=1,NVAR)
      WRITE(9,'(A)')' FIRST, SECOND, AND THIRD-ORDER GRADIENTS, ETC'
      WRITE(9,*)(GRAD(I),I=1,NVAR)
      WRITE(9,*)(GROLD(I),I=1,NVAR)
      WRITE(9,*)(GROLD2(I),I=1,NVAR)
      WRITE(9,*)(PARREF(I),I=1,NVAR)
      I=ILOOP+1
      WRITE(9,*)ETOT,ESCF,EKIN,DELOLD,DELTAT,DLOLD2,I,
     +TOTIME,OLDT,TREF,OLDH,OLDT,REFSCF
         LINEAR=(NORBS*(NORBS+1))/2
         WRITE(10)(PA(I),I=1,LINEAR)
         IF(NALPHA.NE.0)WRITE(10)(PB(I),I=1,LINEAR)
      WRITE(6,'(//10X,'' RUNNING OUT OF TIME, RESTART FILE WRITTEN'')')
      STOP
      ENDIF
  230 CONTINUE
      END
