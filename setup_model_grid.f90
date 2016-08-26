  SUBROUTINE setup_model_grid() 

  ! Set up model grid cells

  USE types

  IMPLICIT NONE    

  INTEGER   :: I, J, K, L

  IF (model_type .EQ. 1) THEN
     ! Read 1D wind model
     CALL read_1D_model()
  ELSE IF (model_type .EQ. 2) THEN
     CALL read_2D_model()
  ELSE IF (model_type .EQ. 3) THEN
     ! Read 3D wind model
     !CALL read_3D_model()
     PRINT*, 'ERROR: Unknown model type', model_type
     STOP
  ELSE
     PRINT*, 'ERROR: Unknown model type', model_type
     STOP
  END IF

  END SUBROUTINE setup_model_grid
