SUBROUTINE read_3D_nico()

USE types
IMPLICIT NONE

INTEGER                                 :: cur_mgi, n_mgi
DOUBLE PRECISION                        :: cur_Teff, cur_rinf, cur_vinf
DOUBLE PRECISION                        :: x, y, z, vx, vy, vz, rho, temp
DOUBLE PRECISION                        :: cur_radius, cur_velocity, cur_rstar
DOUBLE PRECISION                        :: cur_xmax, cur_ymax, cur_zmax

INTEGER                                 :: reading_models

CHARACTER(500)                           :: line, modelfile

INTEGER                                  :: numbions, atom_number
INTEGER                                  :: I, J



modelfile=TRIM(inputmodelFile)

add_mg = 1
cur_mgi = 0

OPEN(UNIT=11, FILE=modelfile)
 n_mgi = 0
 ! number of cells
 DO
  READ(11, '(A)', IOSTAT=reading_models) line
  IF(line(1:1) == '*') CYCLE
  IF(reading_models /= 0) EXIT

  n_mgi = n_mgi + 1

 END DO

 n_modelgrid = n_mgi

 ALLOCATE(model_grid(n_modelgrid + add_mg))

 REWIND(11)

 cur_Teff = 1.e20
 cur_rstar = 1.D40
 cur_rinf = 1.e0
 cur_vinf = 1.e0
 cur_xmax = 0.e0
 cur_ymax = 0.e0
 cur_zmax = 0.e0

 ! reading data
 DO

  READ(11, '(A)', IOSTAT=reading_models) line
  write(*,*) 'read_3D_nico: line = ', line
  IF(line(1:1) == '*') CYCLE
  IF(reading_models /= 0) EXIT

  cur_mgi = cur_mgi + 1

  READ(line,*) x, y, z, vx, vy, vz, rho, temp
  
  model_grid(cur_mgi)%vec_pos(1) = x * R_sun
  model_grid(cur_mgi)%vec_pos(2) = y * R_sun
  model_grid(cur_mgi)%vec_pos(3) = z * R_sun

  cur_radius = R_sun * sqrt(x**2 + y**2 + z**2)

  model_grid(cur_mgi)%vec_vel(1) = vx
  model_grid(cur_mgi)%vec_vel(2) = vy
  model_grid(cur_mgi)%vec_vel(3) = vz

  cur_velocity = sqrt(vx**2 + vy**2 + vz**2)

  model_grid(cur_mgi)%rho = rho
  model_grid(cur_mgi)%T = temp

  IF(temp < cur_Teff) cur_Teff = temp
  IF(cur_radius < cur_rstar) cur_rstar = cur_radius
  IF(cur_radius > cur_rinf) cur_rinf = cur_radius
  IF(cur_velocity > cur_vinf) cur_vinf = cur_velocity

  IF(abs(x) > cur_xmax) cur_xmax = abs(x)
  IF(abs(y) > cur_ymax) cur_ymax = abs(y)
  IF(abs(z) > cur_zmax) cur_zmax = abs(z)

  
  model_grid(cur_mgi)%J = 0.D0
  model_grid(cur_mgi)%assoc_cells = 0

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

 T_eff = cur_Teff
 R_star = cur_rstar
 R_inf = cur_rinf
 V_inf = cur_vinf
 xmax = cur_xmax * R_sun
 ymax = cur_ymax * R_sun
 zmax = cur_zmax * R_sun

CLOSE(11)


END SUBROUTINE read_3D_nico
