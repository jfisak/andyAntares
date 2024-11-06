SUBROUTINE connect_3D_pseudo()

USE types
IMPLICIT NONE

INTEGER                                 :: cur_propcell





DO cur_propcell = 1, n_propgcells
 dyn_cell(cur_propcell)%model_index = cur_propcell
 model_grid(cur_propcell)%assoc_cells = model_grid(cur_propcell)%assoc_cells + 1
END DO














END SUBROUTINE connect_3D_pseudo
