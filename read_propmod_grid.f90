! reads saved propmod grid from a file
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE read_propmod_grid()

USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: ind_I, ind_J
CHARACTER(LEN=filename_length)                  :: propmod_file

DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos, vel, width
DOUBLE PRECISION                                :: rho, temp, volume, el_dens

INTEGER                                         :: down_cell, up_cell, model_index
INTEGER, DIMENSION(6)                           :: neighbors

INTEGER                                         :: Nx, Ny, Nz

INTEGER                                         :: atom_number, numbions, assoc_cells

LOGICAL                                         :: is_diff

write(*,*) '**************************************'
write(*,*) 'READING THE PROPMOD GRID FROM THE FILE'
write(*,*) '**************************************'

! now it works only for 3D models!!!
write(propmod_file,"(A, A12)") TRIM(outputfolder), '/propmod.dat'

! OPEN(49, FORM='unformatted', FILE=propmod_file)
OPEN(49, FILE=propmod_file)

! saves the model grid
read(49, *) n_modelgrid
read(49, *) T_eff
read(49, *) R_star
read(49, *) R_inf
read(49, *) V_inf
read(49, *) xmax, ymax, zmax
read(49, *) Nx, Ny, Nz
read(49, *) add_mg

write(*,*) 'read_propmod_grid: n_modelgrid = ', n_modelgrid
ALLOCATE(model_grid(n_modelgrid + add_mg))

nx_cell = Nx
ny_cell = Ny
nz_cell = Nz
Ngrid = nx_cell * ny_cell * nz_cell

DO ind_I = 1, n_modelgrid
 read(49, *) pos, vel, rho, temp, volume, assoc_cells, el_dens, is_diff
 write(*,*) 'read_propmod_grid: pos = ', pos

 IF(model_type == 1) THEN
  model_grid(ind_I)%rwind = pos(ind_x)
  model_grid(ind_I)%vel = vel(ind_x)
 ELSE IF(model_type == 3) THEN
  model_grid(ind_I)%vec_pos = pos
  model_grid(ind_I)%rwind = pos(ind_x)
  model_grid(ind_I)%angle = pos(ind_y)
  model_grid(ind_I)%vec_vel = vel
 END IF
 model_grid(ind_I)%rho = rho
 model_grid(ind_I)%T = temp
 model_grid(ind_I)%volume = volume
 model_grid(ind_I)%assoc_cells = assoc_cells
 model_grid(ind_I)%J = 0.D0
 model_grid(ind_I)%e_dens = el_dens
 model_grid(ind_I)%is_difapp = is_diff

 ALLOCATE (model_grid(ind_I)%grid_comp(n_elements))
 DO ind_J = 1, n_elements
  numbions = elements(ind_J)%nions
  ALLOCATE (model_grid(ind_I)%grid_comp(ind_J)%grid_ion(numbions))
  atom_number = elements(ind_J)%atom_number
  model_grid(ind_I)%grid_comp(ind_J)%abund = elements(ind_J)%abundance
 END DO
END DO

model_grid(n_modelgrid+1)%rwind = 0.D0
model_grid(n_modelgrid+1)%vel   = 0.D0
model_grid(n_modelgrid+1)%rho   = 0.D0


read(49, *) n_propgcells
write(*,*) 'read_propmod_grid. n_propgcells = ', n_propgcells
ALLOCATE(dyn_cell(n_propgcells))
DO ind_I = 1, n_propgcells
 read(49, *) pos, width, up_cell, down_cell, neighbors, model_index
 
 dyn_cell(ind_I)%corner = pos
 dyn_cell(ind_I)%width = width
 dyn_cell(ind_I)%up_cell = up_cell
 dyn_cell(ind_I)%down_cell = down_cell
 dyn_cell(ind_I)%neighbor = neighbors
 dyn_cell(ind_I)%model_index = model_index
END DO

basic_cell_width = dyn_cell(1)%width


CLOSE(49)






















END SUBROUTINE read_propmod_grid
