      SUBROUTINE DEFINE_EVENT_FLAG_CLUSTER

C PROCEDURE TO CREATE A COMMON EVENT FLAG CLUSTER AVAILABLE TO THE CURRENT
C PROCESS AND OTHERS COMMUNICATING WITH IT
C CHECKS TO SEE IF FLAG 64 IS SET (BIT 64), PROCESS CHECKS THIS FLAG
C AND CAN BE DIRECTED TO A NEW PATH

C *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *   *
C       WRITTEN BY DARYL S. BOUDREAUX 5-JUL-82                         *
C             ENERGY AND SPECIALTIES MATERIALS                          *
C             CORPORATE TECHNOLOGY                                      *
C             ALLIED CORPORATION                                        *
C             COLUMBIA ROAD                                             *
C             MORRISTOWNSHIP, N.J. 07960                                *
C                                                                       *
C      REVISED SLIGHTLY ON 16-JUL-82 BY J. E. WILKALIS                  *
C      VMS 3.0 UPDATE BY CAPT DM STORCH, 04 NOV 82                      *
C                                                                       *
C  NO WARRANTY, EXPRESS OR IMPLIED, IS MADE BY THE AUTHOR OR ALLIED     *
C  CORPORATION AS TO THE ACCURACY OR FUNCTIONING OF THE PROGRAM, AND NO *
C  RESPONSIBILITY IS ASSUMED BY THESE PARTIES IN CONNECTION THEREWITH!  *
C                                                                       *
C *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *   *

      INTEGER STATUS,             ! RETURN CODE FOR SYSTEM SERVICE CALLS
     1	ID_FLAG 
      CHARACTER JOB_NAME*8,NAME_SUBSTRING*8,NAME*15

C FOLLOWING COMMON IS USED FOR $GETJPI SYSTEM ROUTINE
      COMMON/ITEMLIST/
     1 W_LEN1,W_CODE1,L_ADDR1,L_LENADDR1,
     2 W_LEN2,W_CODE2,L_ADDR2,L_LENADDR2,
     3 W_LEN3,W_CODE3                    !W_CODE3=0 SIGNIFIES THE END

	CHARACTER*1 NULL
	DATA NULL/0/

C      PARAMETER      ! DEFINES SOME SYSTEM SERVICE RETURN CODES
C     +    SS$_NORMAL='00000001'X,
C     +    SS$_WASSET='00000009'X,
C     +    SS$_WASCLR='00000001'X

	EXTERNAL SS$_NORMAL,SS$_WASSET,SS$_WASCLR

C DEFINE NAMES OF SYSTEM SERVICE FUNCTIONS:
      INTEGER SYS$ASCEFC          ! CREATES COMMON EVENT FLAG CLUSTER

C OBTAIN PROCESS ID AND JOB NAME USING SYSTEM ROUTINE $GETJPI
      CALL PROCESS_INFO(JOB_PID,JOB_NAME)
C# 1000 FORMAT(1X,'DEFINED EVENT FLAG CLUSTER FOR JOB ',A )

	LL= INDEX(JOB_NAME,NULL) - 1
C#	WRITE (6,*) 'LENGTH OF NAME=',LL
      LENGTH=LEN(JOB_NAME)
	IF (LENGTH .GT. LL) LENGTH = LL
	IF (LENGTH .LE. 0) LENGTH = 8
C#      WRITE(6,1000)JOB_NAME(:LENGTH)
      IF (LENGTH.GT.8) LENGTH=8 ! UP TO 8 CHARACTERS OF SPECIFIED NAME ARE USED
      NAME_SUBSTRING=JOB_NAME(1:LENGTH)
      NAME=NAME_SUBSTRING//'_EFLAGS' ! IS FULL NAME OF CLUSTER


C     ASSUME THAT ONLY ONE LONGWORD OF FLAGS (32) IS NEEDED AND USE THOSE
C     IN CLUSTER 2, NUMBERED 64-95 (LOWER NUMBERS ARE LOCAL EVENT FLAGS)

      STATUS=SYS$ASCEFC(%VAL(64),NAME,,)
      IF ( STATUS .NE. %LOC(SS$_NORMAL) ) THEN
        WRITE(6,100) NAME
  100 FORMAT(' FAILURE IN ATTEMPT TO CREATE COMMON EVENT FLAG CLUSTER: ',
     +A15)
      STOP
	ELSE
	   ID_FLAG=65
	   STATUS = SYS$SETEF(%VAL(ID_FLAG))
	   IF (.NOT. STATUS ) RETURN
	   WRITE (6,*) 'ERROR IN SETTING ID FLAG, STATUS=',STATUS
	ENDIF
      RETURN
      END
 
      SUBROUTINE CHECK_EVENT_FLAG(NUMBER,SET)

C     IF EVENT FLAG: NUMBER IS SET, RETURNS SET=.TRUE. OTHERWISE SET=.FLASE.

      LOGICAL SET
      INTEGER STATUS,STATE ! RETURN CODE FOR SYSTEM SERVICE CALLS AND STATE
                           ! OF EVENT FLAG CLUSTER WHEN EXAMINED

C      PARAMETER      ! DEFINES SOME SYSTEM SERVICE RETURN CODES
C     +    SS$_NORMAL='00000001'X,
C     +    SS$_WASSET='00000009'X,
C     +    SS$_WASCLR='00000001'X

	EXTERNAL SS$_NORMAL,SS$_WASSET,SS$_WASCLR

C     DEFINE NAMES OF SYSTEM SERVICE FUNCTIONS:
      INTEGER SYS$READEF! READ COMMON EVENT FLAG CLUSTER
C
      STATUS=SYS$READEF(%VAL(NUMBER),STATE)
      IF    (STATUS.EQ. %LOC(SS$_WASSET) ) THEN 
                                     SET=.TRUE.
                                     RETURN
      ELSEIF(STATUS.EQ. %LOC(SS$_WASCLR) ) THEN
                                     SET=.FALSE.
                                     RETURN
      ELSE 
        WRITE(6,100) NUMBER
  100   FORMAT(' UNSUCCESSFUL ATTEMPT TO READ EVENT FLAG NUMBER',I3)
        STOP
      ENDIF
      END


      SUBROUTINE PROCESS_INFO(JOB_PID,JOB_NAME)

C THIS ROUTINE USES THE SYSTEM SERVICE ROUTINE $GETJPI TO OBTAIN
C MANY TYPES OF INFO ON THE PROCESS (IN THIS CASE SEE PARAM LIST)

C SET UP DATA TYPES W=WORD, L=LONGWORD

      IMPLICIT INTEGER*2(W)
      IMPLICIT INTEGER*4(L)
      INTEGER*4 STATUS
      CHARACTER*8 JOB_NAME

C CREATE A FILE WITH THE FOLLOWING MACRO PROGRAM
C      $JPIDEF GLOBAL
C      .END
C ASSEMBLE IT USING:
C $MACRO/LIST FILE_NAME
C THIS GIVES YOU A LISTING OF THE PARAMETER VALUES REQUIRED TO OBTAIN
C THE INFORMATION DESIRED

C      PARAMETER JPI$_PID    = '00000319'X
C      PARAMETER JPI$_PRCNAM = '0000031C'X

	EXTERNAL JPI$_PID,JPI$_PRCNAM

C DECLARE WORKING STORAGE
      INTEGER*4 BUFFER_VALUES(3)
      INTEGER*4 LENGTH_VALUES(2)

C DECLARE SYS$GETJPI ITEM LIST AS A COMMON
      COMMON/ITEMLIST/
     1 W_LEN1,W_CODE1,L_ADDR1,L_LENADDR1,
     2 W_LEN2,W_CODE2,L_ADDR2,L_LENADDR2,
     3 W_LEN3,W_CODE3                    !W_CODE3=0 SIGNIFIES END

C INITIALIZE ALL STATIC VALUES IN THE ITEM LIST
      DATA W_LEN1/4/,W_LEN2/8/,W_LEN3/4/
C      DATA W_CODE1/JPI$_PID/
C      DATA W_CODE2/JPI$_PRCNAM/
      DATA W_CODE3/0/                    !END OF INFO DESIRED

	W_CODE1 = %LOC(JPI$_PID)
	W_CODE2 = %LOC(JPI$_PRCNAM)

C ITEM FIELDS REQUIRING ADDRESSES MUST BE ASSIGNED AT RUN TIME
      L_ADDR1=%LOC(BUFFER_VALUES(1))
      L_ADDR2=%LOC(BUFFER_VALUES(2))
      L_LENADDR1=%LOC(LENGTH_VALUES(1))
      L_LENADDR2=%LOC(LENGTH_VALUES(2))

C PERFORM THE SYSTEM SERVICE CALL
      CALL SYS$GETJPI(,,,W_LEN1,STATUS,,)
C      STATUS = SYS$GETJPI(,,,W_LEN1,,,)!THIS CALL DID NOT WORK

      IF(STATUS.EQ.1)THEN
          JOB_PID=BUFFER_VALUES(1)
          J=2 ! IS LONGWORD COUNTER IN BUFFER_VALUES
          DO I=1,LENGTH_VALUES(2)
             IF (I.EQ.5) J=3 ! 5'TH BYTE IN NEXT LONGWORD
             JOB_NAME(I:I)=CHAR(BUFFER_VALUES(J))

C SHIFT BUFFER AFTER EACH LOOP TO OBTAIN NEXT CHARACTER
             BUFFER_VALUES(J)=BUFFER_VALUES(J)/256
          END DO
      ELSE
         WRITE(6,10)STATUS
      ENDIF
   10 FORMAT(' PROCESS_INFO NOT OBTAINED',I8)
      RETURN
      END
