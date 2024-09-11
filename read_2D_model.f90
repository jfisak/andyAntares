! subroutine for reading a 2D model grid
SUBROUTINE read_2D_model()

USE types
USE constants

IMPLICIT NONE
 
 SELECT CASE (inputModel)
  ! reading data from Petr Kurfurst disc model
  ! these files are in this form:
  ! 1. radius / m
  ! 2. coordinate in perpendicular direction to the
  !     disc plane / m
  ! 3. density / kg/m^3
  ! 4. radial velocity / m/s
  ! 5. angular velocity / m/s
  ! 6. temperature / K
  CASE(1)
   CALL read_2D_basic()
  CASE(2)
   CALL read_2D_peku()
  CASE DEFAULT
   STOP 'unknown type of 2D model'
 END SELECT

END SUBROUTINE read_2D_model
