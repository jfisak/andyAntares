SUBROUTINE save_output()

USE types
IMPLICIT NONE


! folder variables
CHARACTER(LEN=30)                       :: outputfolder
LOGICAL                                 :: dirExists
CHARACTER(LEN=30)                       :: mkdirCMD
CHARACTER(LEN=60)                       :: lineOutput
! save informations about lines
INTEGER                                 :: I
DOUBLE PRECISION                        :: wavle

! creates a folder, where an output will be saved
! it reads a shell variable OUTPUTFO, if it does not
! exist, it will create (or not, if it already exists)
! a directory 3dwindmodel
CALL GET_ENVIRONMENT_VARIABLE("OUTPUTFO", outputfolder)
IF(outputfolder(:) == '') THEN
 outputfolder = '3dwindmodel'
END IF

inquire( file=trim(outputfolder)//'/.', exist=dirExists )

IF(.NOT. dirExists) THEN
 mkdirCMD = 'mkdir '//TRIM(outputfolder)
 CALL SYSTEM(mkdirCMD)
 write(*,*) 'save_output: creating a folder: ', outputfolder
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now will save important variables
! line rates
lineOutput = trim(outputfolder)//'/linevar.dat'
write(*,*) 'save_output: lineOutput = ', TRIM(lineOutput)
OPEN(11,FILE=lineOutput)
 DO I = 1, ntransitions
  ! wavelength is in Angstroms
  wavle = 1e8 * light_speed / linelist(I)%freq
  WRITE(11,*) elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, wavle, linelist(I)%n_int
 END DO
CLOSE(11)





END SUBROUTINE save_output
