      SUBROUTINE SURFCITY(COORD,JFLAG)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      CHARACTER KEYWRD*80
      LOGICAL JFLAG
      COMMON /OPTIM / IMP,IMP0,LEC,IPRT
      COMMON /KEYWRD/ KEYWRD
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /BORN  / BP(NUMATM),FGB(NPACK),CCT1,ZEFF(NUMATM),
     1                QEFF(NUMATM)
      COMMON /SURF  / SURFCT,SURFACT(NUMATM),ATAR(NUMATM),ITYPE(NUMATM)
      COMMON /SMOOTH/ RAL(NUMATM)
      dimension COORD(3,*),RA(NUMATM),BORN(NUMATM),WAALS(100)
      dimension t_cos(90),t_sin(90),surfen1(100),chion(100)
      dimension xp(4000),yp(4000),zp(4000),surfen2(100)
      dimension angs(45),hits(4000),qflec(100)
      dimension distsq(5000),R(NPACK),wools(100),rmax(numatm)
c
      data radmax/-100.0/
c
C     VAN DER WAALS RADII A LA BONDI, JPC 68 (1964) 441.
C
      data waals/
     1 1.20,                                                      0.00,
     2 1.82,0.00,                       0.00,1.70,1.55,1.52,1.47 ,1.54,
     3 2.27,1.73,                       2.50,2.10,1.80,1.80,1.75 ,1.88,
     4 2.75,0.00,7*0.00,1.63,1.40,1.39, 2.40,2.10,1.85,1.90,1.85 ,2.02,
     5 0.00,0.00,7*0.00,1.63,1.72,1.58, 2.50,2.20,2.10,2.06,1.98 ,2.60,
     6  46*0.00/
c
C     COULOMB INTEGRAL BASES (PSEUDO RADII) FOR GENERATION OF BORN
C     RADII AND INTERACTIONS. (see below)
C
      data wools/
     1 0.57,                                                      0.00,
     2 1.82,0.00,                       0.00,1.68,1.40,1.46,1.37 ,1.54,
     3 2.27,1.73,                       2.50,2.10,1.80,1.30,1.65 ,1.88,
     4 2.75,0.00,7*0.00,1.63,1.40,1.39, 2.40,2.10,1.85,1.90,1.75 ,2.02,
     5 0.00,0.00,7*0.00,1.63,1.72,1.58, 2.50,2.20,2.10,2.06,1.88 ,2.60,
     6  46*0.00/
C
C     CHARGE ACCOUNTING ARCTANGENT PRE-FACTOR. THUS, EFFECTIVE COULOMB
C     RADIUS = WOOLS + CHION * f(ARCTAN(f(CHARGE)))
C
      data chion/
     1   1.303,                                                   0.00,
     2   0.00 , 0.00,          0.00,  0.00,  0.62, -0.25,  0.181, 0.00,
     3   0.00 , 0.00,          0.00,  0.00,  0.00,  0.80,  0.618, 0.00,
     4   0.00 , 0.00, 10*0.00, 0.00,  0.00,  0.00,  0.00,  0.705, 0.00,
     5   0.00 , 0.00, 10*0.00, 0.00,  0.00,  0.00,  0.00,  0.932, 0.00,
     6   46*0.00/
C
C    CHARGE REPRESENTING THE INFLECTION POINT ON THE COULOMB INTEGRAL
C    CURVE DESCRIBED AS AN ARCTANGENT FUNCTION
C
      data qflec/
     1   -0.3 ,                                                   0.00,
     2   0.00 , 0.00,          0.00,  0.50,  0.40,  0.75,  0.700, 0.00,
     3   0.00 , 0.00,          0.00,  0.00,  0.00,  0.70,  0.750, 0.00,
     4   0.00 , 0.00, 10*0.00, 0.00,  0.00,  0.00,  0.00,  0.700, 0.00,
     5   0.00 , 0.00, 10*0.00, 0.00,  0.00,  0.00,  0.00,  0.600, 0.00,
     6   46*0.00/
C
C     SURFACE AREA DEPENDENT CORRECTION TERMS DETERMINED WITH NO
C     ACCOUNT TAKEN OF CHEMICAL ENVIRONMENT (AQUO)
C
      data surfen1/
     1   0.00,                                                    0.00,
     2   0.00,  0.00,          0.00, 14.95,-73.65,-52.74, 18.47,  0.00,
     3   0.00,  0.00,          0.00,  0.00,  0.00,-18.80, -2.14,  0.00,
     4   0.00,  0.00, 10*0.00, 0.00,  0.00,  0.00,  0.00, -9.11 , 0.00,
     5   0.00,  0.00, 10*0.00, 0.00,  0.00,  0.00,  0.00, -8.21 , 0.00,
     6   46*0.00/
C
C     SURFACE AREA DEPENDENT CORRECTION TERMS DETERMINED WITH
C     ACCOUNT TAKEN OF CHEMICAL ENVIRONMENT (ENVAQ)
C
C     CURRENTLY: C(H), -nH, -oH, -sH,
C              sp3/Amide-N, sp2/Arom/sp-N, sp3-O, sp2-O,
C              F, S, Cl, Br, I (total of 13 types)
C
      data surfen2/  4.15, 58.96, -23.39, 49.49, -368.97,
     1 -47.38,-109.70,-25.61, 21.17,-44.25,-2.84,-8.93,-13.42,
     2  87*0.00/
      SAVE
      SOLVRAD=1.4D0
      PI=DACOS(-1.D0)
      CELEID=1.D0/78.3D0
      IFLAG=0
      IEND=0
      TSTEP=0.05D0
      IF(JFLAG) THEN
            NPTS=50
            NARCS=25
      ELSE
            npts=90
            narcs=45
      ENDIF
       IF(JFLAG) THEN
C     write(6,'(" NUMAT = ",i3)')numat
      do 901 i=1,NUMAT
      SSX=-(DATAN(1.D01*(QEFF(I)+QFLEC(NAT(I))))-PI/2.D0)/PI
      RA(I)=WOOLS(NAT(I))+CHION(NAT(I))*SSX
      IF(RAL(I).EQ.0) THEN
        RAL(I)=RA(I)
      ELSE
        RA(I)=(RA(I)+RAL(I))/2.D0
        RAL(I)=RA(I)
      ENDIF
C     if(.not.jflag) go to 901
C     write(6,'(" radius of atom ",i3," is ",f10.6)')i,ra(i)
901   BORN(I)=0.D0
      ELSE
      do 977 i=1,numat
977   ra(i)=waals(nat(i))+solvrad
      IF(INDEX(KEYWRD,'AQUO').NE.0) THEN
            DO 902 I=1,NUMAT
902         SURFACT(I)=SURFEN1(NAT(I))
      ELSE
            DO 903 I=1,NUMAT
903         SURFACT(I)=SURFEN2(ITYPE(I))
      ENDIF
       ENDIF
c
C     CALCULATE INTERATOMIC DISTANCES
C
905   DO 15 I=1,NUMAT
      RMAX(I)=0.d0
      DO 15 J=1,I
        IJ=(I*(I-1))/2+J
        XX1=COORD(1,I)-COORD(1,J)
        XX2=COORD(2,I)-COORD(2,J)
        XX3=COORD(3,I)-COORD(3,J)
      R(IJ)=SQRT(XX1*XX1+XX2*XX2+XX3*XX3)
15    rmax(i)=max(rmax(i),r(ij))
c
c           Calculate the constants needed to define and rotate
c           the points
c
            nrcs=narcs
            rotang=(PI/2.D0)/((nrcs+1)/2.D0)
            nacd2=nrcs/2
            do 20 i=(-nacd2),(nacd2)
20          angs(i+nacd2+1)=rotang*i
C
C     THIS MASTER LOOP CALCULATES TOTAL EXPOSED AREA OF SPHERE I.
C
      do 80 iat=1,nUMat
      IF(.NOT.JFLAG) GO TO 910
C
C     SAVE ORIGINAL RA IN AHOLD AND ADJUST FIRST SPHERE TO COULOMB SURFACE
C
          AHOLD=RA(IAT)
 700      XAREA=4.0d0*PI*RA(IAT)*RA(IAT)
c
c           Zero the summation variables
c
910         ipts=0
            ipsum=0
c
c           Now, evaluate the points on the surface of the sphere
c
            do 50 ia=1,nrcs
c
c                 Construct the dots on the surface of the sphere
c
                  zdisp=ra(iat)*sin(angs(ia))
                  rad  =ra(iat)*cos(angs(ia))
c
c                 Calculate the number of points in the circle
c
                        np=int((rad/ra(iat))*npts)
                        xterm=(2.D0*PI)/np
c
c                       Build the list of theta's
c
                        do 30 j=1,np
                              theta=(j-1)*xterm
                              t_cos(j)=cos(theta)
30                      t_sin(j)=sin(theta)
c
c                       Update the list of points
c
                        do 40 j=ipsum+1,ipsum+np
                         xp(j)=t_cos(j-ipsum)*rad+COORD(1,iat)
                         yp(j)=t_sin(j-ipsum)*rad+COORD(2,iat)
                         zp(j)=zdisp             +COORD(3,iat)
40                      hits(j)=0.0
                  ipsum=ipsum+np
50    CONTINUE
c
c           Point vector is ready - process it
c
            itz=iat-1
            do 65 j=1,itz
            wid=ra(iat)+ra(j)
            lm=(iat*(iat-1))/2+j
            if(wid.lt.r(lm)) go to 65
                  do 60 i=1,ipsum
                     distsq(i)=((COORD(1,j)-xp(i))**2+
     1                     (COORD(2,j)-yp(i))**2+
     2                     (COORD(3,j)-zp(i))**2)-
     3                     ra(j)**2
60                hits(i)=min(hits(i),distsq(i))
65          continue
            izt=iat+1
            do 68 j=izt,numat
            wid=ra(iat)+ra(j)
            lm=(j*(j-1))/2+iat
            if(wid.lt.r(lm)) go to 68
                  do 67 i=1,ipsum
                     distsq(i)=((COORD(1,j)-xp(i))**2+
     1                     (COORD(2,j)-yp(i))**2+
     2                     (COORD(3,j)-zp(i))**2)-
     3                     ra(j)**2
67                hits(i)=min(hits(i),distsq(i))
68          continue
c
c           Count the visible points
c
            do 70 i=1,ipsum
70          ipts=ipts+dsign(1.d0,1.d0*hits(i))+1
            ipts=ipts/2
            area=(4.D0*PI*ra(iat)*RA(IAT)*IPTS/IPSUM)
            if(.not.jflag) then
                  atar(iat)=area
                  surfact(iat)=surfact(iat)*area/1.d3
            else
                        PRCNT=AREA/XAREA
C
C     IFLAG DETERMINES WHETHER THE CALCULATION IS FOR THE CENTER OF A
C     BORN SHELL OR THE INNER SURFACE. THE FINAL SHELL MUST ENTIRELY 
C     ENCLOSE THE REMAINDER OF THE MOLECULE OR SUPERMOLECULE. IEND FLAGS
C     THIS CONDITION HAVING BEEN MET
C
                  IF(IFLAG.EQ.1) GO TO 704
                        IF(RA(IAT).LT.RMAX(IAT)) GO TO 703
                        IF(PRCNT.LT.0.999) GO TO 703
                        IEND=1
                  IF(TSTEP.EQ.0.05D0) TSTEP=0.D0
703               IFLAG=1
                        RA(IAT)=RA(IAT)+TSTEP
                        GO TO 700
704               IFLAG=0
                        IF(IEND.EQ.1) GO TO 705
            TDIFF=1.d0/(RA(IAT)-TSTEP)-1.d0/(RA(IAT)+TSTEP)
            BORN(IAT)=BORN(IAT)+PRCNT*TDIFF
                        RA(IAT)=RA(IAT)+TSTEP
                        TSTEP=1.5d0*TSTEP
                        GO TO 700
C
C     ADD THE FINAL SHELL'S CONTRIBUTION AND RESET ALL THE STEP VARIABLES
C     TO INITIAL VALUES. RESET THE PROPER RADIUS FOR ATOM I.
C
705               BORN(IAT)=BORN(IAT)+1.d0/RA(IAT)
                        IEND=0
                        TSTEP=0.05d0
                        RA(IAT)=AHOLD
      endif
80    CONTINUE
      IF(.NOT.JFLAG) THEN
                  surfct=0.d0
            DO 710 M=1,NUMAT
710         surfct=surfct+surfact(m)
      ELSE
                  DO 714 M=1,NUMAT
714         BORN(M)=1.d0/BORN(M)
c
C     COMPUTE THE STILL FACTOR FGB AS A FUNCTION OF R AND ALPHA
C
                  DO 888 I=1,NUMAT
                  DO 888 J=1,I
                  IJ=(I*(I-1))/2+J
                  K=MIN(NAT(I),NAT(J))
                  L=MAX(NAT(I),NAT(J))
                  KL=(L*(L-1))/2+K
                  OOFAC=0.D0
                  IF(KL.EQ.22) THEN
                        IF(R(IJ).LE.2.399.AND.R(IJ).GE.1.401) THEN
                              COGH=4.d0
                              COGW=-1.75D0
                              TRIP=((R(IJ)-1.9D0)/0.5D0)**2
                              OOFAC=COGH*EXP(COGW/(1-TRIP))
                              GO TO 831
                        ENDIF
                  ENDIF
                  IF(KL.EQ.36) THEN
                        IF(R(IJ).LE.2.599.AND.R(IJ).GE.0.801) THEN
                              COGH=9.d0
                              COGW=-1.75D0
                              TRIP=((R(IJ)-1.7D0)/0.9D0)**2
                              OOFAC=COGH*EXP(COGW/(1-TRIP))
                        ENDIF
                  ENDIF
831               AIJ=BORN(I)*BORN(J)
                  DD=EXP(-R(IJ)*R(IJ)/(4.D0*AIJ))+OOFAC
                        FGB(IJ)=SQRT(R(IJ)*R(IJ)+AIJ*DD)
888         CONTINUE
      ENDIF
      RETURN
      END
