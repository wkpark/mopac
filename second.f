      FUNCTION SECOND()
      DOUBLE PRECISION SECOND
C******************************************************
C
C   SECOND, ON EXIT, CONTAINS THE NUMBER OF CPU SECONDS
C   SINCE THE START OF THE CALCULATION.
C
C******************************************************
      LOGICAL SETOK
      CHARACTER*1 X
      DIMENSION A(2)
      DATA SETOK   /  .TRUE.    /, SHUT/0.D0/
      Y=ETIME(A)
      CPU=A(1)
***********************************************************************
*
*   NOW TO SEE IF A FILE LOGICALLY CALLED SHUTDOWN EXISTS, IF IT DOES
*   THEN INCREMENT CPU TIME BY 1,000,000 SECONDS.
*
************************************************************************
      OPEN(UNIT=4, FILE='SHUTDOWN',STATUS='UNKNOWN')
      READ(4,'(A)',END=10, ERR=10)X
*
*          FILE EXISTS, THEREFORE INCREMENT TIME
*
      SHUT=1.D6
      IF( SETOK) THEN
         WRITE(6,'(///10X,''****   JOB STOPPED BY OPERATOR   ****'')')
         SETOK=.FALSE.
      ENDIF
   10 CONTINUE
      SECOND=CPU+SHUT
      RETURN
      END
