! reads electron density from a file
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE read_e_nd()

USE types
USE constants

IMPLICIT NONE

DOUBLE PRECISION                          :: junk
CHARACTER(filename_length)                             :: modelfile

DOUBLE PRECISION                                :: eldens
INTEGER                                         :: ind_I

write(*,*) '******************************************'
write(*,*) 'READING THE ELECTRON DENSITY FROM THE FILE'
write(*,*) '******************************************'

modelfile=TRIM(inputmodelFile)

IF(model_type == 1) THEN
 SELECT CASE (inputModel)
  CASE(0)
   OPEN(UNIT=11, FILE=modelfile)
    READ(11,*) junk
    READ(11,*) junk
    READ(11,*) junk
    READ(11,*) junk
    DO ind_I = 1, n_modelgrid
     READ(11,*) junk, junk, junk, junk, junk, eldens
     model_grid(ind_I)%e_dens = eldens
     ! write(*,*) 'read_e_nd: ind_I = ', I, ' eldens = ', eldens
    END DO
   CLOSE(11)
 END SELECT
ELSE IF(model_type == 3) THEN
 SELECT CASE(inputModel)
 CASE(1)
  
 END SELECT
END IF

END SUBROUTINE read_e_nd
