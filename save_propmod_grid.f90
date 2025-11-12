! this sbr saves propagation and model grid into a single file for a purpose to read
! it later instead of creating it in every run
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE save_propmod_grid()

USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: ind_I
CHARACTER(LEN=filename_length)                  :: propmod_file

INTEGER                                         :: n_propgrid

DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos, vel, width
DOUBLE PRECISION                                :: rho, temp, volume, el_dens

INTEGER                                         :: down_cell, up_cell, model_index, assoc_cells
INTEGER, DIMENSION(6)                           :: neighbors

LOGICAL                                         :: is_diff


write(propmod_file,"(A, A12)") TRIM(outputfolder), '/propmod.dat'

OPEN(49, FILE=propmod_file)

! saves the model grid
write(49,*) n_modelgrid
write(49,*) T_eff
write(49,*) R_star
write(49,*) R_inf
write(49,*) V_inf
write(49,*) xmax, ymax, zmax
write(49,*) nx_cell, ny_cell, nz_cell
write(49,*) add_mg

DO ind_I = 1, n_modelgrid
 IF(model_type == 1) THEN
  pos(ind_x) = model_grid(ind_I)%rwind
  pos(ind_y) = 0.e0
  pos(ind_z) = 0.e0
  vel(ind_x) = model_grid(ind_I)%vel
  vel(ind_y) = 0.e0
  vel(ind_z) = 0.e0
 ELSE IF(model_type == 2) THEN
  pos(ind_x) = model_grid(ind_I)%rwind
  pos(ind_y) = model_grid(ind_I)%angle
  pos(ind_z) = 0.e0
  vel(ind_x) = model_grid(ind_I)%vel
  vel(ind_y) = model_grid(ind_I)%velang
  vel(ind_z) = 0.e0
 ELSE IF(model_type == 3) THEN
  pos = model_grid(ind_I)%vec_pos
  vel = model_grid(ind_I)%vec_vel
 END IF
 rho = model_grid(ind_I)%rho
 temp = model_grid(ind_I)%T
 volume = model_grid(ind_I)%volume
 assoc_cells = model_grid(ind_I)%assoc_cells
 el_dens = model_grid(ind_I)%e_dens
 is_diff = model_grid(ind_I)%is_difapp

 write(49,*) pos, vel, rho, temp, volume, assoc_cells, el_dens, is_diff
END DO

! saves the associated propagation grid
write(49,*) n_propgcells
DO ind_I = 1, n_propgcells
 pos = dyn_cell(ind_I)%corner
 width = dyn_cell(ind_I)%width
 up_cell = dyn_cell(ind_I)%up_cell
 down_cell = dyn_cell(ind_I)%down_cell
 neighbors = dyn_cell(ind_I)%neighbor
 model_index = dyn_cell(ind_I)%model_index

 write(49,*) pos, width, up_cell, down_cell, neighbors, model_index
END DO










CLOSE(49)


















END SUBROUTINE save_propmod_grid

