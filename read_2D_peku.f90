! Petr Kurfurst Supernova model
! column; physical quantity; unities
! 1: radial coordinate, m
! 2: angle coordinate, rad
! 3: mass density, kg/m^3
! 4: radial velocity, m/s
! 5: polar velocity, rad/s
! 6: temperature, K
SUBROUTINE read_2D_peku()

USE types
USE constants

IMPLICIT NONE


DOUBLE PRECISION                          :: junk
INTEGER, PARAMETER                     :: maxrows = 60000000

DOUBLE PRECISION                        :: radius, temp, dens, velrad, velang, angle
INTEGER                                 :: atom_number, numbions
INTEGER                                 :: ios
INTEGER                                 :: I, J

CHARACTER(len=400)                      :: ch_line

DOUBLE PRECISION, PARAMETER             :: min_temp=5.E4


add_mg = 2
write(99,*) 'we will read input input data from Petr Kurfurst model of stellar disc'
! firstly we calculate number of rows in the file
n_modelgrid = 0
OPEN(UNIT=15,status='old', FILE=inputmodelFile)
 READ(15, *) junk
 READ(15, *) junk
 DO I = 1, maxrows
  READ(15,'(A)',IOSTAT = ios) ch_line
  if(ch_line == '') cycle
  if(ios /= 0) EXIT
  ! testing of the temperature: points with too low temperature are excluded
  READ(ch_line,*) radius, angle, dens, velrad, velang, temp
  if(temp < min_temp) cycle
 n_modelgrid = n_modelgrid + 1
 END DO
 write(99,*) 'mumber of model grids: ', n_modelgrid
 IF (n_modelgrid .EQ. 0) STOP 'no model grid cells were found...'
 ! n_modelgrid + 1 ... for dummy cells
 ! n_modelgrid + 2 ... for cells with r < R_inf but too far from some model grid point
 !                     (vacuum cell) 
 ALLOCATE ( model_grid(n_modelgrid + add_mg))
 REWIND(15)
 READ(15, *) T_eff
 READ(15, *) R_star
 I = 1
 DO 
  READ(15,'(A)',IOSTAT = ios) ch_line
  if(ios /= 0) EXIT
  if(ch_line == '') cycle
  ! write(*,*) 'read_2d_model: ch_line = ', ch_line
  READ(ch_line,*) radius, angle, dens, velrad, velang, temp
  if(temp < min_temp) cycle
  model_grid(I)%rwind = radius * 1.D2
  model_grid(I)%angle = angle
  model_grid(I)%vel = velrad * 1.D2
  model_grid(I)%velang = velang * 1.D2
  model_grid(I)%rho = dens * 1.D-3
  model_grid(I)%T = temp
  model_grid(I)%J = 0.D0
  model_grid(I)%assoc_cells = 0
  ! write(*,*) 'read_2d_model: I = ', I, ' temp = ', temp
  ALLOCATE(model_grid(I)%grid_comp(n_elements))
  ! now we add informations about every included element for every model cell
  DO J = 1, n_elements
   numbions = elements(J)%nions
   ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
   atom_number = elements(J)%atom_number
   model_grid(I)%grid_comp(J)%abund = elements(J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
   !model_grid(I)%grid_comp(J)%numb_den = tot_nd
  END DO
  I = I + 1
 END DO
 R_inf = MAXVAL(model_grid(:)%rwind)
 ! write(*,*) 'read_2d_model: R_inf = ', R_inf, ' R_star = ', R_star
CLOSE(15)



END SUBROUTINE read_2D_peku
