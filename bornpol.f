      SUBROUTINE BORNPOL(NITER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON /OPTIM / IMP,IMP0,LEC,IPRT
      COMMON /DENSTY/ P(MPACK)
      COMMON /GEOM  / GEO(3,NUMATM)
      COMMON /MOLKST/ NUMAT,NAT(NUMATM),NFIRST(NUMATM),NMIDLE(NUMATM),
     1                NLAST(NUMATM), NORBS, NELECS,NALPHA,NBETA,
     2                NCLOSE,NOPEN,NDUMY,FRACT
      COMMON /BORN  / BP(NUMATM),FGB(NPACK),CCT1,ZEFF(NUMATM),
     1                QEFF(NUMATM)
      COMMON /SURF  / SURFCT,SURFACT(NUMATM),ATAR(NUMATM),ITYPE(NUMATM)
      DIMENSION COORD(3,NUMATM),ZHELP(100)
      DATA ZHELP/
     1 1.,                                                         2.,
     2 1., 2.,                                 3., 4., 5., 6., 7., 8.,
     3 1., 2.,                                 3., 4., 5., 6., 7., 8.,
     4 1., 2.,        10*2.,                   3., 4., 5., 6., 7., 8.,
     5 1., 2.,        10*2.,                   3., 4., 5., 6., 7., 8.,
     6 1., 2.,        24*2.,                   3., 4., 5., 6., 7., 8.,
     7 1., 2.,        12*2./
        SAVE
C
C       THIS SUBROUTINE CALCULATES THE BORN POLARIZATION
C       CONTRIBUTION TO THE (I,I)'TH FOCK MATRIX ELEMENT
C       FGB IS A PACKED MATRIX OF STILL FGB FACTORS
C       SEE SUBROUTINE SURFCITY FOR REFERENCES
C
C       THE MATRIX BP CONTAINS THE ATOM SPECIFIC ADJUSTMENTS TO THE DIAGONAL
C       ELEMENTS OF THE FOCK MATRIX AND IS DETERMINED SELF-CONSISTENTLY
C
C       THE CHARGE CORRECTION TERM CCT1 REFLECTS THE INTERACTION OF THE
C       "EFFECTIVE" NUCLEAR CHARGES WITH THE INDUCED SOLVATON CHARGES
C
      CELEID=1.d0/78.3d0
      NATOMS=NUMAT
      CCT1=0.d0
      K=0
      IF(ZEFF(1).GT.0) GO TO 81
C
C       CALCULATE THE "EFFECTIVE" NUCLEAR CHARGE--I.E. VALENCE PROTONS
C
      DO 55 I=1,NATOMS
55    ZEFF(I)=ZHELP(NAT(I))
C
C       CALCULATE ATOMIC CHARGES
C
81    DO 57 I=1,NATOMS
        BP(I)=0.D0
        QEFF(I)=ZEFF(I)
        DO 57 J=NFIRST(I),NLAST(I)
        K=K+J
57      QEFF(I)=QEFF(I)-P(K)
      IF(SQRT(FLOAT(NITER)).NE.INT(SQRT(FLOAT(NITER)))) GO TO 73
C     write(6,'(" Iteration ",i3)')niter
      CALL GMETRY(GEO,COORD)
      CALL SURFCITY(COORD,.TRUE.)
      CALL SURFCITY(COORD,.FALSE.)
C        call vecprt(fgb,numat)
73    DO 100 I=1,NATOMS
        JP=I+1
        DO 120 J=1,I
        IJ=(I*(I-1))/2+J
        BP(I)=BP(I)+QEFF(J)/FGB(IJ)
120     CCT1=CCT1-7.1983d0*(1.d0-CELEID)*ZEFF(J)*QEFF(I)/FGB(IJ)
        DO 130 K=JP,NATOMS
        IJ=(K*(K-1))/2+I
        BP(I)=BP(I)+QEFF(K)/FGB(IJ)
130     CCT1=CCT1-7.1983d0*(1.d0-CELEID)*ZEFF(K)*QEFF(I)/FGB(IJ)
100   continue
      DO 200 I=1,NATOMS
        BP(I)=14.3966d0*(1.d0-CELEID)*BP(I)
200   continue
      CCT1=CCT1+SURFCT/23.06d0
      RETURN
      END
