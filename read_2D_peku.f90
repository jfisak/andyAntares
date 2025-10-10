! Petr Kurfurst Supernova model
! column; physical quantity; unities
! 1: radial coordinate, r
! 2: vertical coordinate, z
! 3: mass density, kg/m^3
! 4: radial velocity, m/s
! 5: polar velocity, rad/s
! 6: temperature, K
SUBROUTINE read_2D_peku()

USE types
USE constants

IMPLICIT NONE


CHARACTER(150)                          :: junk
INTEGER, PARAMETER                     :: maxrows = 60000000

DOUBLE PRECISION                        :: radius, temp, dens, velrad, velz, coor_z
INTEGER                                 :: atom_number, numbions
INTEGER                                 :: ios
INTEGER                                 :: ind_I, ind_J

CHARACTER(len=filename_length)                      :: ch_line

DOUBLE PRECISION, PARAMETER             :: min_temp=1.5E4


add_mg = 3
write(99,*) 'we will read input input data from Petr Kurfurst model of stellar disc'
! firstly we calculate number of rows in the file
n_modelgrid = 0
OPEN(UNIT=15,status='old', FILE=inputmodelFile)
 READ(15, *) junk
 READ(15, *) junk
 DO ind_I = 1, maxrows
  READ(15,'(A)',IOSTAT = ios) ch_line
  if(ch_line == '') cycle
  if(ios /= 0) EXIT
  ! testing of the temperature: points with too low temperature are excluded
  READ(ch_line,*) radius, coor_z, dens, temp, junk, velrad, velz
  if(temp < min_temp) cycle
  n_modelgrid = n_modelgrid + 1
 END DO
 write(99,*) 'mumber of model grids: ', n_modelgrid
 ! write(*,*) 'mumber of model grids: ', n_modelgrid
 IF (n_modelgrid .EQ. 0) STOP 'no model grid cells were found...'
 ! n_modelgrid + 1 ... for dummy cells
 ! n_modelgrid + 2 ... for cells with r < R_inf but too far from some model grid point
 !                     (vacuum cell) 
 outerspace_index = n_modelgrid + 1
 vacuum_index = n_modelgrid + 2
 photosphere_index = n_modelgrid + 3

 ALLOCATE ( model_grid(n_modelgrid + add_mg))
 REWIND(15)
 READ(15, *) T_eff
 READ(15, *) junk
 ind_I = 1
 DO 
  READ(15,'(A)',IOSTAT = ios) ch_line
  if(ios /= 0) EXIT
  if(ch_line == '') cycle
  ! write(*,*) 'read_2d_model: ch_line = ', ch_line
  READ(ch_line,*) radius, coor_z, dens, temp, junk, velrad, velz
  if(temp < min_temp) cycle
  model_grid(ind_I)%rxywind = radius
  model_grid(ind_I)%zwind = coor_z
  model_grid(ind_I)%rwind = sqrt(radius**2+coor_z**2)
  model_grid(ind_I)%vel = velrad
  model_grid(ind_I)%velz = velz
  model_grid(ind_I)%rho = dens
  model_grid(ind_I)%T = temp
  model_grid(ind_I)%J = 0.D0
  model_grid(ind_I)%assoc_cells = 0
  if(isnan(velrad)) STOP 'read_2d_model: velrad = NaN'
  ALLOCATE(model_grid(ind_I)%grid_comp(n_elements))
  ! now we add informations about every included element for every model cell
  DO ind_J = 1, n_elements
   numbions = elements(ind_J)%nions
   ALLOCATE (model_grid(ind_I)%grid_comp(ind_J)%grid_ion(numbions))
   atom_number = elements(ind_J)%atom_number
   model_grid(ind_I)%grid_comp(ind_J)%abund = elements(ind_J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
   !model_grid(I)%grid_comp(J)%numb_den = tot_nd
  END DO
  ind_I = ind_I + 1
 END DO
 R_star = 1.D2
 R_inf = MAXVAL(model_grid(:)%rwind)
 model_grid(outerspace_index)%vel = 0.D0
 model_grid(vacuum_index)%vel = 0.D0
 model_grid(outerspace_index)%velz = 0.D0
 model_grid(vacuum_index)%velz = 0.D0
CLOSE(15)



END SUBROUTINE read_2D_peku
