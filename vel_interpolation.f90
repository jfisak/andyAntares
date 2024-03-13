SUBROUTINE vel_interpolation()

USE types
IMPLICIT NONE

INTEGER                                         :: n_prop_grid, cur_propgrid_cell


! going through every propGrid cell
DO cur_propgrid_cell = 1, n_propgcells
 IF(dyn_cell(cur_propgrid_cell)%up_cell == 0) THEN
  
 END IF ! dyn_cell(cur_propgrid_cell)%up_cell == 0
END DO






END SUBROUTINE vel_interpolation
