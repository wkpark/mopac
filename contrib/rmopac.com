$!   COMMAND FILE TO RUN A MOPAC JOB. THIS FILE SHOULD RESIDE IN THE
$!   MOPAC DIRECTORY. IT SHOULD BE ACCESSED FROM MOPAC.COM, BUT CAN
$!   STAND ALONE, IF NECESSARY.
$!
$!   THE CALL IS	$RMOPAC filename  directory  
$!
$	SHOW TIME
$!
$	SET VERIFY
$!
$!	Make assignments
$!
$DEL*ETE :== DELETE
$	IF P2 .NES. ""	THEN SET DEFAULT 'P2'
$	OPEN /ERR=NOEND1 DUMMY 'P1'.END  ! see if shutdown file is there
$       DELETE/NOCONFIRM 'P1'.END;*
$NOEND1:
$	ASSIGN 'P1'.DAT FOR005
$	ASSIGN 'P1'.OUT FOR006
$	ASSIGN 'P1'.RES FOR009
$	ASSIGN 'P1'.DEN FOR010
$	ASSIGN 'P1'.ARC FOR012
$	ASSIGN 'P1'.GPT FOR013
$       ASSIGN 'P1'.END SHUTDOWN
$!
$	ON ERROR	THEN GOTO AA
$	ON CONTROL_Y	THEN GOTO AA	! Cleanup if ^Y
$	RUN MOPACDIRECTORY:MOPAC
$!	SET NOVERIFY
$!
$	SHOW TIME
$!
$!	Delete assignments
$!
$ AA:
$	SET NOCONTROL_Y			! Continue cleanup if ^Y
$       OPEN /ERR=NOEND DUMMY 'P1'.END 	! see if shutdown file is there
$       DELETE/NOCONFIRM 'P1'.END;*
$NOEND: CLOSE DUMMY
$	SET CONTROL_Y
$	SET NOVERIFY
$!
$!	END
$!
