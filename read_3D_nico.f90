SUBROUTINE read_3D_nico()

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: cur_mgi, n_mgi
DOUBLE PRECISION                        :: cur_rinf, cur_vinf
DOUBLE PRECISION                        :: x, y, z, vx, vy, vz, rho, temp, lambda
DOUBLE PRECISION                        :: cur_radius, cur_velocity
DOUBLE PRECISION                        :: cur_xmax, cur_ymax, cur_zmax

INTEGER                                 :: reading_models

CHARACTER(500)                           :: line, modelfile

INTEGER                                  :: numbions, atom_number
INTEGER                                  :: I, J



modelfile=TRIM(inputmodelFile)

! add_mg = 1 r < R_star
! add_mg = 2 r > R_inf
! add_mg = 3 r > R_star && r < R_inf, vacuum cell
add_mg = 3
cur_mgi = 0

OPEN(UNIT=11, FILE=modelfile)
 n_mgi = 0
 ! number of cells
 READ(11, '(A)', IOSTAT=reading_models) line
 READ(11, '(A)', IOSTAT=reading_models) line
 READ(11, '(A)', IOSTAT=reading_models) line
 DO
  READ(11, '(A)', IOSTAT=reading_models) line
  IF(line(1:1) == '*') CYCLE
  IF(reading_models /= 0) EXIT

  n_mgi = n_mgi + 1

 END DO

 n_modelgrid = n_mgi
 write(*,*) 'read_3D_nico: n_modelgrid = ', n_modelgrid, ' add_mg = ', add_mg

 ALLOCATE(model_grid(n_modelgrid + add_mg))

 REWIND(11)

 READ(11,*) T_eff
 READ(11,*) R_star
 READ(11,*) R_inf
 cur_rinf = 1.e0
 cur_vinf = 1.e0
 cur_xmax = 0.e0
 cur_ymax = 0.e0
 cur_zmax = 0.e0

 ! reading data
 DO

  READ(11, '(A)', IOSTAT=reading_models) line
  IF(line(1:1) == '*') CYCLE
  IF(reading_models /= 0) EXIT

  cur_mgi = cur_mgi + 1

  READ(line,*) x, y, z, vx, vy, vz, rho, temp, lambda
  
  model_grid(cur_mgi)%vec_pos(1) = x * R_sun
  model_grid(cur_mgi)%vec_pos(2) = y * R_sun
  model_grid(cur_mgi)%vec_pos(3) = z * R_sun
  model_grid(cur_mgi)%rwind = sqrt(x**2 + y**2 + z**2) * R_sun

  cur_radius = R_sun * sqrt(x**2 + y**2 + z**2)

  model_grid(cur_mgi)%vec_vel(1) = vx
  model_grid(cur_mgi)%vec_vel(2) = vy
  model_grid(cur_mgi)%vec_vel(3) = vz

  model_grid(cur_mgi)%diff_param = lambda

  cur_velocity = sqrt(vx**2 + vy**2 + vz**2)
  IF(cur_velocity > light_speed) THEN
   write(*,*) 'read_3D_nico: the velocity of the point I = ', cur_velocity,&
    ' is larger than the speed of light'
   CALL abort()
  END IF

  model_grid(cur_mgi)%rho = rho
  model_grid(cur_mgi)%T = temp

  IF(cur_radius > cur_rinf) THEN
   ! cur_rinf = cur_radius
   cur_vinf = cur_velocity
  END IF

  IF(abs(x) > cur_xmax) cur_xmax = abs(x)
  IF(abs(y) > cur_ymax) cur_ymax = abs(y)
  IF(abs(z) > cur_zmax) cur_zmax = abs(z)

  
  model_grid(cur_mgi)%J = 0.D0
  model_grid(cur_mgi)%assoc_cells = 0
  model_grid(cur_mgi)%volume = 0.D0

  ALLOCATE (model_grid(cur_mgi)%grid_comp(n_elements))
  DO J = 1, n_elements
   numbions = elements(J)%nions
   ALLOCATE (model_grid(cur_mgi)%grid_comp(J)%grid_ion(numbions))
   atom_number = elements(J)%atom_number
   !model_grid(I)%grid_comp(J)%abund = massfrac(atom_number)
   model_grid(cur_mgi)%grid_comp(J)%abund = elements(J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass
   !model_grid(I)%grid_comp(J)%numb_den = tot_nd
  END DO

 END DO

! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
! All cells out of model grid set to 0 and associate to n_modelgrid.
! Other cells will obtain particular values with memory
DO I = 1, add_mg
 model_grid(n_modelgrid + I)%vec_vel = (/ 0.e0, 0.e0, 0.e0 /)
 model_grid(n_modelgrid + I)%rwind = 0.D0
 model_grid(n_modelgrid + I)%vel   = 0.D0
 model_grid(n_modelgrid + I)%rho   = 0.D0
 model_grid(n_modelgrid + I)%J = 0.D0
 model_grid(n_modelgrid + I)%assoc_cells = 0
 model_grid(n_modelgrid + I)%volume = 0.D0
END DO
  
 V_inf = cur_vinf/10.0
 xmax = (cur_xmax + cur_xmax / R_star / 100.0) * R_sun  
 ymax = (cur_ymax + cur_ymax / R_star / 100.0) * R_sun 
 zmax = (cur_zmax + cur_zmax / R_star / 100.0) * R_sun 
 ! write(*,*) 'read_3D_nico: T_eff = ', T_eff, 'R_star = ', R_star/R_inf, ' R_inf = ', R_inf/R_sun, ' V_inf = ', V_inf
 ! write(*,*) 'read_3D_nico: xmax = ', xmax/R_inf, ' ymax = ', ymax/R_inf, ' zmax = ', zmax/R_inf
 !  write(*,*) 'read_3D_nico: xmax/R_inf = ', xmax/R_inf
 write(*,*) 'read_3D_nico: V_prop/V_inf = ', 6.0*(xmax/R_inf)*(ymax/R_inf)*(zmax/R_inf)/const_pi

CLOSE(11)


END SUBROUTINE read_3D_nico
