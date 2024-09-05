SUBROUTINE read_propmod_grid()

USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: I, J
CHARACTER(LEN=60)                               :: propmod_file

INTEGER                                         :: n_propgrid

DOUBLE PRECISION, DIMENSION(3)                  :: pos, vel, width
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

DO I = 1, n_modelgrid
 read(49, *) pos, vel, rho, temp, volume, assoc_cells, el_dens, is_diff

 IF(model_type == 1) THEN
  model_grid(I)%rwind = pos(1)
  model_grid(I)%vel = vel(1)
 ELSE IF(model_type == 3) THEN
  model_grid(I)%vec_pos = pos
  model_grid(I)%rwind = norm2(pos)
  model_grid(I)%vec_vel = vel
 END IF
 model_grid(I)%rho = rho
 model_grid(I)%T = temp
 model_grid(I)%volume = volume
 model_grid(I)%assoc_cells = assoc_cells
 model_grid(I)%J = 0.D0
 model_grid(I)%e_dens = el_dens
 model_grid(I)%is_difapp = is_diff

 ALLOCATE (model_grid(I)%grid_comp(n_elements))
 DO J = 1, n_elements
  numbions = elements(J)%nions
  ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
  atom_number = elements(J)%atom_number
  model_grid(I)%grid_comp(J)%abund = elements(J)%abundance
 END DO
END DO

model_grid(n_modelgrid+1)%rwind = 0.D0
model_grid(n_modelgrid+1)%vel   = 0.D0
model_grid(n_modelgrid+1)%rho   = 0.D0


read(49, *) n_propgcells
write(*,*) 'read_propmod_grid. n_propgcells = ', n_propgcells
ALLOCATE(dyn_cell(n_propgcells))
DO I = 1, n_propgcells
 read(49, *) pos, width, up_cell, down_cell, neighbors, model_index
 
 dyn_cell(I)%corner = pos
 dyn_cell(I)%width = width
 dyn_cell(I)%up_cell = up_cell
 dyn_cell(I)%down_cell = down_cell
 dyn_cell(I)%neighbor = neighbors
 dyn_cell(I)%model_index = model_index
END DO

basic_cell_width = dyn_cell(1)%width


CLOSE(49)






















END SUBROUTINE read_propmod_grid
