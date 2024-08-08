SUBROUTINE vel_pseudo3D_model()

USE types
IMPLICIT NONE

INTEGER                                         :: cur_pgi, cur_mgi
DOUBLE PRECISION, DIMENSION(3)                  :: cur_centre

DO cur_pgi = 1, n_propgcells
 cur_mgi = dyn_cell(cur_pgi)%model_index
 IF(cur_mgi <= n_modelgrid) THEN
  dyn_cell(cur_pgi)%vec_vel = model_grid(cur_mgi)%vec_vel
 ELSE IF(cur_mgi == n_modelgrid + 1) THEN
  dyn_cell(cur_pgi)%vec_vel = (/ 0.D0, 0.D0, 0.D0 /)
 ELSE IF(cur_mgi == n_modelgrid + 2) THEN
  cur_centre = dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width/2.0
  dyn_cell(cur_pgi)%vec_vel = R_inf * cur_centre/norm2(cur_centre)
 END IF

END DO










END SUBROUTINE vel_pseudo3D_model
