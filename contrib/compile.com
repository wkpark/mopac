$ASS DIMSIZES.DAT SIZES
$!*******************************************************************
$!                                                                  *
$! COM FILE TO COMPILE AND LINK MOPAC GIVEN A FILE CALLED MOPAC.OPT *
$!                                                                  *
$! ADAPTED BY JJPS FROM A COM FILE WRITTEN BY DONN STORCH, USAFA    *
$!*******************************************************************
$!
$!   THE FORTRAN FILES ARE NOW TO BE COMPILED. 
$A1:
$ AUTO = ""
$ METHOD:== FORTRAN
$ IF P1.EQS."" THEN GOTO COMP
$ IF P1.EQS."LINK" THEN GOTO A
$ LINE = P1
$ GOTO GET_MODULE
$COMP:
$ ON ERROR THEN GOTO COMP_ABORT
$!
$!  AUTOMATIC COMPILATION VIA MOPAC.OPT
$!
$ OPEN/READ FREDD MOPAC.OPT
$READ_LINE:
$ READ/ERR=FINISHED/END=FINISHED FREDD LINE
$GET_MODULE:
$ IF F$LEN(LINE).LT.1 THEN GOTO READ_LINE
$ IF F$EXT(0,1,LINE).NES." " THEN GOTO GOT_SOMETHING
$ LINE = F$EXT(1,F$LEN(LINE)-1,LINE)
$ GOTO GET_MODULE
$GOT_SOMETHING:
$ IF F$EXT(0,1,LINE).EQS."-" THEN GOTO READ_LINE
$ IF F$EXT(0,1,LINE).EQS."DZRO_M" THEN GOTO READ_LINE
$ L2=F$LOC(",",LINE)
$ P1 = F$EXT(0,L2,LINE)
$ IF F$LEN(AUTO).EQ.0 THEN GOTO COMPILE_IT
$ IF P1.NES.AUTO THEN GOTO SKIP
$ AUTO = ""
$COMPILE_IT:
$ WRITE SYS$OUTPUT "''METHOD' ''P1'"
$ 'METHOD' 'P1'
$SKIP:
$ LINE = F$EXT(L2+1,F$LEN(LINE),LINE)
$ GOTO GET_MODULE
$FINISHED:
$ IF P1.EQS."" THEN CLOSE FREDD
$ GOTO A
$!
$! HERE THE OBJECT FILES ARE ALL LINKED TOGETHER.
$!
$A:
$! INQUIRE P1 "Do you want to link MOPAC [Y/N] (N)
$! IF P1.NES."Y" GOTO B
$ LINK MOPAC/OPT
$!
$B:
$ INQUIRE P1 "DO YOU WANT TO CREATE THE MOPAC HELP LIBRARY (Y/N) "
$	IF P1.NES."Y" THEN GOTO D
$LIBRARY/CREATE/HELP MOPAC MOPAC
$! here I discover the next available library name
$ IF "''F$LOG( "HLP$LIBRARY")'" .NES. "" THEN GOTO START_LOOP
$ WRITE SYS$OUTPUT "Logical name HLP$LIBRARY should be used for MOPAC.HLB"
$ GOTO D
$START_LOOP:
$ LIB_NR = 0
$LIB_LOOP:
$ LIB_NR = 'LIB_NR' + 1
$ TTEMPP:== HLP$LIBRARY_'LIB_NR'
$ IF "''F$LOG( TTEMPP )'" .NES. "" THEN GOTO LIB_LOOP
$ WRITE SYS$OUTPUT "Logical name ''TTEMPP' should be used for MOPAC.HLB"
$D:
$ WRITE SYS$OUTPUT ""
$ IF "''F$LOG( "MOPACHELP" )'" .NES. "" THEN GOTO E
$ WRITE SYS$OUTPUT "Be sure to define the logical symbol MOPACHELP."
$E:
$ WRITE SYS$OUTPUT ""
$	WRITE SYS$OUTPUT -
	" COMPILATION COMMAND FILE FINISHED "
$ SET NOVERIFY
$ EXIT
$COMP_ABORT:
$ SET NOVERIFY
$ EXIT
