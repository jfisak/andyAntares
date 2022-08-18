SUBROUTINE read_3D_pseudo3D()


USE types
IMPLICIT NONE

CHARACTER(80)                                   :: modelfile
CHARACTER(80)                                   :: junk

INTEGER                                         :: reading_grid

DOUBLE PRECISION                                :: dens, temp
DOUBLE PRECISION, DIMENSION(3)                  :: pos, vel

INTEGER                                         :: I, J
INTEGER                                         :: numbions, atom_number
INTEGER                                         :: Nx, Ny, Nz

modelfile=TRIM(inputmodelFile)

OPEN(UNIT=11, FILE=modelfile)
 READ(11,*) T_eff
 READ(11,*) R_star
 READ(11,*) R_inf
 READ(11,*) V_inf
 READ(11,*) xmax, ymax, zmax
 READ(11,*) Nx, Ny, Nz

 add_mg = 1

 xmax = (xmax + xmax / Nx) * R_star
 ymax = (ymax + ymax / Ny) * R_star
 zmax = (zmax + zmax / Nz) * R_star

 write(99,*) 'setting new propagation grid size to: nx_cell =', nx_cell, ' ny_cell = ', ny_cell, &
  ' nz_cell = ', nz_cell
 nx_cell = Nx
 ny_cell = Ny
 nz_cell = Nz

 DO
  READ(11, *, iostat = reading_grid) junk
  IF(reading_grid /= 0) EXIT
  n_modelgrid = n_modelgrid + 1
 END DO
 ALLOCATE (model_grid(n_modelgrid + add_mg))
 REWIND(11)
 READ(11,*) junk
 READ(11,*) junk
 READ(11,*) junk
 READ(11,*) junk
 READ(11,*) junk
 READ(11,*) junk

 DO I = 1, n_modelgrid
  READ(11, *) pos, vel, dens, temp
  model_grid(I)%vec_pos = pos  * R_star
  model_grid(I)%vec_vel = vel
  model_grid(I)%rho = dens
  model_grid(I)%T = temp
  model_grid(I)%J = 0.D0
  model_grid(I)%assoc_cells = 0
  
  ALLOCATE (model_grid(I)%grid_comp(n_elements))
  DO J = 1, n_elements
   numbions = elements(J)%nions
   ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
   atom_number = elements(J)%atom_number
   !model_grid(I)%grid_comp(J)%abund = massfrac(atom_number)
   model_grid(I)%grid_comp(J)%abund = elements(J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass
   !model_grid(I)%grid_comp(J)%numb_den = tot_nd
  END DO

 END DO

  ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
  ! All cells out of model grid set to 0 and associate to n_modelgrid.
  ! Other cells will obtainde particular values with memory
  model_grid(n_modelgrid+1)%rwind = 0.D0
  model_grid(n_modelgrid+1)%vel   = 0.D0
  model_grid(n_modelgrid+1)%rho   = 0.D0



CLOSE(11)

END SUBROUTINE read_3D_pseudo3D
