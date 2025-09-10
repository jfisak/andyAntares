! reads testing pseudo 3D model
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE read_3D_pseudo3D()


USE types
USE constants
IMPLICIT NONE

CHARACTER(filename_length)                                   :: modelfile
CHARACTER(filename_length)                                   :: junk

INTEGER                                         :: reading_grid

DOUBLE PRECISION                                :: dens, temp
! DOUBLE PRECISION                                :: mod_xmax, mod_ymax, mod_zmax
! DOUBLE PRECISION                                :: mod_xmin, mod_ymin, mod_zmin
! DOUBLE PRECISION                                :: len_x, len_y, len_z
! DOUBLE PRECISION                                :: width_x, width_y, width_z
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos, vel

INTEGER                                         :: ind_I, ind_J
INTEGER                                         :: numbions, atom_number
INTEGER                                         :: Nx, Ny, Nz
INTEGER, DIMENSION(1)                           :: ind_min
! INTEGER(KIND=8)                                 :: int_xmax, int_ymax, int_zmax
! INTEGER(KIND=8)                                 :: int_xmin, int_ymin, int_zmin


modelfile=TRIM(inputmodelFile)

vel_propgrid = .true.
vel_modgrid = .false.

OPEN(UNIT=11, FILE=modelfile)
 READ(11,*) T_eff
 READ(11,*) R_star
 READ(11,*) R_inf
 READ(11,*) V_inf
 READ(11,*) Nx, Ny, Nz

 add_mg = 3

 write(*,*) 'read_3D_pseudo3D: my_rank = ', my_rank
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
 DO ind_I = 1, n_modelgrid
  READ(11, *) pos, vel, dens, temp
  model_grid(ind_I)%vec_pos = pos * R_star
  model_grid(ind_I)%rwind = norm2(pos) * R_star
  model_grid(ind_I)%vec_vel = vel
  model_grid(ind_I)%vel = norm2(vel)
  model_grid(ind_I)%rho = dens
  model_grid(ind_I)%T = temp
  model_grid(ind_I)%J = 0.D0
  model_grid(ind_I)%assoc_cells = 0
  
  ALLOCATE (model_grid(ind_I)%grid_comp(n_elements))
  DO ind_J = 1, n_elements
   numbions = elements(ind_J)%nions
   ALLOCATE (model_grid(ind_I)%grid_comp(ind_J)%grid_ion(numbions))
   atom_number = elements(ind_J)%atom_number
   !model_grid(ind_I)%grid_comp(ind_J)%abund = massfrac(atom_number)
   model_grid(ind_I)%grid_comp(ind_J)%abund = elements(ind_J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(ind_I)%grid_comp(ind_J)%abund / elements(ind_J)%atom_mass
   !model_grid(ind_I)%grid_comp(ind_J)%numb_den = tot_nd
  END DO

 END DO


 ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
 ! All cells out of model grid set to 0 and associate to n_modelgrid.
 ! Other cells will obtainde particular values with memory
 photosphere_index = n_modelgrid + 1
 outerspace_index = n_modelgrid + 2
 vacuum_index = n_modelgrid + 3
 model_grid(n_modelgrid + add_mg)%vec_vel = (/ 0.e0, 0.e0, 0.e0 /)
 ! cells below the lower boundary
 model_grid(photosphere_index)%rwind = 0.D0
 model_grid(photosphere_index)%vel   = 0.D0
 model_grid(photosphere_index)%rho   = 0.D0
 ! cells beyond the outer boundary
 model_grid(outerspace_index)%rwind = R_inf
 model_grid(outerspace_index)%vel   = V_inf
 model_grid(outerspace_index)%rho   = 0.D0
 ! vacuum cells
 model_grid(vacuum_index)%rwind = 0.D0
 model_grid(vacuum_index)%vel   = 0.D0
 model_grid(vacuum_index)%rho   = 0.D0

 xmin = MINVAL(model_grid(:)%vec_pos(ind_x))
 ymin = MINVAL(model_grid(:)%vec_pos(ind_y))
 zmin = MINVAL(model_grid(:)%vec_pos(ind_z))
 xmax = MAXVAL(model_grid(:)%vec_pos(ind_x))
 ymax = MAXVAL(model_grid(:)%vec_pos(ind_y))
 zmax = MAXVAL(model_grid(:)%vec_pos(ind_z))

 ind_min(1) = MINLOC(model_grid(:)%rwind,1, MASK=(model_grid(:)%rwind >= R_star))
 V_star = model_grid(ind_min(1))%vel
! len_x = mod_xmax - mod_xmin
! len_y = mod_ymax - mod_ymin
! len_z = mod_zmax - mod_zmin
!
! width_x = len_x / (Nx - 1)
! width_y = len_y / (Ny - 1)
! width_z = len_z / (Nz - 1)
!
! int_xmin = CEILING(mod_xmin - width_x / 2.0, kind=8)
! int_ymin = CEILING(mod_ymin - width_y / 2.0, kind=8)
! int_zmin = CEILING(mod_zmin - width_z / 2.0, kind=8)
!
! int_xmax = CEILING(mod_xmax + width_x / 2.0, kind=8)
! int_ymax = CEILING(mod_ymax + width_y / 2.0, kind=8)
! int_zmax = CEILING(mod_zmax + width_z / 2.0, kind=8)
!
! xmin = DBLE(int_xmin)
! ymin = DBLE(int_ymin)
! zmin = DBLE(int_zmin)
! xmax = DBLE(int_xmax)
! ymax = DBLE(int_ymax)
! zmax = DBLE(int_zmax)
 
! write(*,*) 'read_3D_pseudo3D: T_eff = ', T_eff, 'R_star = ', R_star, ' R_inf = ', R_inf, ' V_inf = ', V_inf
! write(*,*) 'read_3D_pseudo3D: xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax
! STOP 'read_3D_pseudo3D: testing'

CLOSE(11)




END SUBROUTINE read_3D_pseudo3D
