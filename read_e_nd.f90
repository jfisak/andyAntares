SUBROUTINE read_e_nd()

USE types
USE constants

IMPLICIT NONE

DOUBLE PRECISION                          :: junk
CHARACTER(80)                             :: modelfile

DOUBLE PRECISION                                :: eldens
INTEGER                                         :: I

modelfile=TRIM(inputmodelFile)

IF(model_type == 1) THEN
 SELECT CASE (inputModel)
  CASE(0)
   OPEN(UNIT=11, FILE=modelfile)
    READ(11,*) junk
    READ(11,*) junk
    READ(11,*) junk
    READ(11,*) junk
    DO I = 1, n_modelgrid
     READ(11,*) junk, junk, junk, junk, junk, eldens
     model_grid(I)%e_dens = eldens
     ! write(*,*) 'read_e_nd: I = ', I, ' eldens = ', eldens
    END DO
   CLOSE(11)
 END SELECT
ELSE IF(model_type == 3) THEN
 SELECT CASE(inputModel)
 CASE(1)
  
 END SELECT
END IF

END SUBROUTINE read_e_nd
