! connects propGrid and modGrid cells
! for each propGrid finds a corresponding modGrid cell
! the variable is saved into dyn_cell(:)%model_index (INT)
! and counted in the model_grid(:)%assoc_cells (INT)
! 1 ... n_modelgrid
! or
! photosphere_index
! outerspace_index
! vacuum_index
!
! INPUT: NONE
! OUTPUT: NONE
!
! the included function will be moved to special subroutines as well as the 3D case
SUBROUTINE connection_prop_model_grid()

USE MPI
USE types
USE constants
USE counters

IMPLICIT NONE
! maximal distance between model and propagation grid
! MUST BE LATER CHANGED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Establish a connection between the propagation grid and the
! model grid. This depends on the model grid type (1D, 2D, 3D)
!__________________________________________________________________________________________________
!__________________________________________________________________________________________________
! 1D case
!__________________________________________________________________________________________________
!__________________________________________________________________________________________________
IF (model_type .EQ. 1) THEN
 ! This is the algorithm needed for a 1D model grid
 ! Define which model grid cell coresponds to the propagation grid cell
 CALL connect_1D_basic()
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ! 2D case
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ELSE IF (model_type .EQ. 2) THEN
   SELECT CASE (inputmodel)
   ! basic 2D model
   CASE(1)
    CALL connect_2D_basic()
   ! PeKu model
   CASE(2)
    CALL connect_2D_peku()
   CASE DEFAULT
    STOP
   END SELECT
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ! 3D case
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ELSE IF (model_type == 3) THEN
   SELECT CASE(inputmodel)
   ! pseudo 3D testing model
   CASE(0)
    CALL connect_3D_pseudo()
   ! hydronico model
   CASE(1)
    CALL connect_3D_hydronico()
   CASE DEFAULT
    write(*,*) 'the choice inputmodel = ', inputmodel, ' is not known'
    STOP
   END SELECT
  END IF
! computing volume of model cells
! DO cur_propcell = 1, max_n_dcell
!  IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
!   gridcell = dyn_cell(cur_propcell)%model_index
!   loc_volume = dyn_cell(cur_propcell)%width(1) * dyn_cell(cur_propcell)%width(2) * dyn_cell(cur_propcell)%width(3)
!   model_grid(gridcell)%volume = model_grid(gridcell)%volume + loc_volume
!   ! counters
!   IF(dyn_cell(cur_propcell)%model_index /= 0) THEN
!    count_pg_mcell = count_pg_mcell + 1
!   ELSE IF(dyn_cell(cur_propcell)%model_index == n_modelgrid + 3) THEN
!    count_pg_vacuum = count_pg_vacuum + 1
!   END IF
!  END IF
! END DO

write(99,*) 'propGrid statistics:'
write(99,*) 'count_pg_mcell = ', count_pg_mcell
write(99,*) 'count_pg_vacuum = ', count_pg_vacuum

END SUBROUTINE connection_prop_model_grid
