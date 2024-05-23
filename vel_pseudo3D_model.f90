SUBROUTINE vel_pseudo3D_model()

USE types
IMPLICIT NONE

INTEGER                                         :: cur_pgi, cur_mgi

DO cur_pgi = 1, n_propgcells
 cur_mgi = dyn_cell(cur_pgi)%model_index
 IF(cur_mgi <= n_modelgrid) THEN
  dyn_cell(cur_pgi)%vec_vel = model_grid(cur_mgi)%vec_vel
 ELSE IF(cur_mgi == n_modelgrid + 1) THEN
  dyn_cell(cur_pgi)%vec_vel = (/ 0.D0, 0.D0, 0.D0 /)
 END IF

END DO










END SUBROUTINE vel_pseudo3D_model
