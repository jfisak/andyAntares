! subroutine for reading a 2D model grid
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE read_2D_model()

USE types
USE constants

IMPLICIT NONE
 
 SELECT CASE (inputModel)
  CASE(1)
   CALL read_2D_basic()
  CASE(2)
   CALL read_2D_peku()
  CASE DEFAULT
   STOP 'unknown type of 2D model'
 END SELECT

END SUBROUTINE read_2D_model
