SUBROUTINE read_3D_pseudo3D()


USE types
USE constants
IMPLICIT NONE

CHARACTER(80)                                   :: modelfile
CHARACTER(80)                                   :: junk

INTEGER                                         :: reading_grid

DOUBLE PRECISION                                :: dens, temp
DOUBLE PRECISION                                :: mod_xmax, mod_ymax, mod_zmax
DOUBLE PRECISION                                :: mod_xmin, mod_ymin, mod_zmin
DOUBLE PRECISION                                :: len_x, len_y, len_z
DOUBLE PRECISION                                :: width_x, width_y, width_z
DOUBLE PRECISION, DIMENSION(3)                  :: pos, vel

INTEGER                                         :: I, J
INTEGER                                         :: numbions, atom_number
INTEGER                                         :: Nx, Ny, Nz


modelfile=TRIM(inputmodelFile)

vel_propgrid = .true.
vel_modgrid = .false.

OPEN(UNIT=11, FILE=modelfile)
 READ(11,*) T_eff
 READ(11,*) R_star
 READ(11,*) R_inf
 READ(11,*) V_inf
 READ(11,*) Nx, Ny, Nz

 add_mg = 2

 write(*,*) 'read_3D_pseudo3D: R_star = ', R_star, ' R_inf = ', R_inf, ' R_inf/R_star = ', R_inf/R_star

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

 write(*,*) 'read_3D_pseudo3D: n_modelgrid = ', n_modelgrid
 DO I = 1, n_modelgrid
  READ(11, *) pos, vel, dens, temp
  model_grid(I)%vec_pos = pos * R_star
  model_grid(I)%rwind = norm2(pos) * R_star
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
 model_grid(n_modelgrid + add_mg)%vec_vel = (/ 0.e0, 0.e0, 0.e0 /)
 ! cells below the lower boundary
 model_grid(n_modelgrid+add_mg - 1)%rwind = 0.D0
 model_grid(n_modelgrid+add_mg - 1)%vel   = 0.D0
 model_grid(n_modelgrid+add_mg - 1)%rho   = 0.D0
 ! cells beyond the outer boundary
 model_grid(n_modelgrid+add_mg)%rwind = R_inf
 model_grid(n_modelgrid+add_mg)%vel   = V_inf
 model_grid(n_modelgrid+add_mg)%rho   = 0.D0

 mod_xmin = MINVAL(model_grid(:)%vec_pos(1))
 mod_ymin = MINVAL(model_grid(:)%vec_pos(2))
 mod_zmin = MINVAL(model_grid(:)%vec_pos(3))
 mod_xmax = MAXVAL(model_grid(:)%vec_pos(1))
 mod_ymax = MAXVAL(model_grid(:)%vec_pos(2))
 mod_zmax = MAXVAL(model_grid(:)%vec_pos(3))

 len_x = mod_xmax - mod_xmin
 len_y = mod_ymax - mod_ymin
 len_z = mod_zmax - mod_zmin

 width_x = len_x / (Nx - 1)
 width_y = len_y / (Ny - 1)
 width_z = len_z / (Nz - 1)

 xmin = mod_xmin - width_x / 2.0
 ymin = mod_ymin - width_y / 2.0
 zmin = mod_zmin - width_z / 2.0

 xmax = mod_xmax + width_x / 2.0
 ymax = mod_ymax + width_y / 2.0
 zmax = mod_zmax + width_z / 2.0
 
  write(*,*) 'read_3D_pseudo3D: T_eff = ', T_eff, 'R_star = ', R_star, ' R_inf = ', R_inf, ' V_inf = ', V_inf
  write(*,*) 'read_3D_pseudo3D: xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax


CLOSE(11)

END SUBROUTINE read_3D_pseudo3D
