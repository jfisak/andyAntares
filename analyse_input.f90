! this sbr analyses input from the file input.dat
SUBROUTINE analyse_input()

USE types
USE constants
IMPLICIT NONE

LOGICAL                         :: inp_model_exists, inputcomp_exists



! test of file existence
! model file
INQUIRE(FILE=inputmodelFile, EXIST=inp_model_exists)

IF(.not. inp_model_exists) THEN
 write(*,*) 'input model file ', inputmodelFile, ' does not exist'
 STOP
END IF

! eldens file
! INQUIRE(FILE=eldensfile, EXIST=eldensfile_exists)
! IF(.not. eldensfile_exists) THEN
!  write(*,*) 'electron density file ', eldensfile, ' does not exist'
!  STOP
! END IF

! input composition
INQUIRE(FILE=inputcomposition, EXIST=inputcomp_exists)
IF(.not. inputcomp_exists) THEN
 write(*,*) 'input composition file ', inputcomposition, ' does not exist'
 STOP
END IF


END SUBROUTINE analyse_input
