! this sbr saves propagation and model grid into a single file for a purpose to read
! it later instead of creating it in every run
SUBROUTINE save_propmod_grid()

USE types
IMPLICIT NONE

INTEGER                                         :: I
CHARACTER(LEN=60)                               :: propmod_file

INTEGER                                         :: n_propgrid

DOUBLE PRECISION, DIMENSION(3)                  :: pos, vel, width
DOUBLE PRECISION                                :: rho, temp, volume

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

DO I = 1, n_modelgrid
 IF(model_type == 1) THEN
  pos(1) = model_grid(I)%rwind
  pos(2) = 0.e0
  pos(3) = 0.e0
  vel(1) = model_grid(I)%vel
  vel(2) = 0.e0
  vel(3) = 0.e0
 ELSE IF(model_type == 3) THEN
  pos = model_grid(I)%vec_pos
  vel = model_grid(I)%vec_vel
 END IF
 rho = model_grid(I)%rho
 temp = model_grid(I)%T
 volume = model_grid(I)%volume
 assoc_cells = model_grid(I)%assoc_cells
 is_diff = model_grid(I)%is_difapp

 write(49,*) pos, vel, rho, temp, volume, assoc_cells, is_diff
END DO

! saves the associated propagation grid
n_propgrid = SIZE(dyn_cell)
write(49,*) n_propgrid
DO I = 1, n_propgrid
 pos = dyn_cell(I)%corner
 width = dyn_cell(I)%width
 up_cell = dyn_cell(I)%up_cell
 down_cell = dyn_cell(I)%down_cell
 neighbors = dyn_cell(I)%neighbor
 model_index = dyn_cell(I)%model_index

 write(49,*) pos, width, up_cell, down_cell, neighbors, model_index
END DO










CLOSE(49)


















END SUBROUTINE save_propmod_grid

