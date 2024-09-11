SUBROUTINE read_2D_basic()

USE types
USE constants


IMPLICIT NONE

INTEGER                                :: ios
INTEGER, PARAMETER                     :: maxrows = 60000000
DOUBLE PRECISION                          :: junk
! loop variables
INTEGER                                :: I, J
! readen physical quantities
DOUBLE PRECISION                       :: radius, perpend, dens, velrad, velang, temp
INTEGER                                :: atom_number, numbions
! calcultaion a stellar radius
DOUBLE PRECISION                       :: act_radius, min_radius, max_radius
DOUBLE PRECISION                       :: act_z, max_z
INTEGER                                :: max_radius_index
INTEGER                                :: vacuum

DOUBLE PRECISION                        :: unit_length, unit_velocity, unit_density

unit_length = 12.64759321736591 * R_sun
unit_velocity = 1.D8
unit_density = 1.41314878888978971775872825028550241D-0006

 add_mg = 2
 write(99,*) 'we will read input input data from the basic 2D model'
 ! firstly we calculate number of rows in the file
 n_modelgrid = 0
  T_eff = 30000
 OPEN(UNIT=15,status='old', FILE=inputmodelFile)
  DO I = 1, maxrows
   READ(15,*,IOSTAT = ios) junk, junk, junk, junk, junk, junk
   if(ios /= 0) EXIT
   if(I == maxrows) THEN
    write(99,*) 'maximum number of records exceeded in subroutine read_2d_model'
    write(99,*) 'exiting program now...'
    STOP
   end if
  n_modelgrid = n_modelgrid + 1
 END DO
 write(99,*) 'mumber of model grids: ', n_modelgrid
 IF (n_modelgrid .EQ. 0) STOP 'no model grid cells were found...'
 ! n_modelgrid + 1 ... for dummy cells
 ! n_modelgrid + 2 ... for cells with r < R_inf but too far from some model grid point
 !                     (vacuum cell) 
 ALLOCATE(model_grid(n_modelgrid + add_mg))
 ! will define vacuum index
 vacuum = n_modelgrid + add_mg
 REWIND(15)
 DO I = 1, n_modelgrid
  READ(15,*) radius, perpend, dens, velrad, velang, temp
  model_grid(I)%rwind = radius * unit_length
  model_grid(I)%angle = perpend
  model_grid(I)%vel = velrad * unit_velocity
  model_grid(I)%velang = velang * unit_velocity
  model_grid(I)%rho = dens * unit_density
  model_grid(I)%T = temp
  model_grid(I)%J = 0.D0
  model_grid(I)%assoc_cells = 0
  ALLOCATE (model_grid(I)%grid_comp(n_elements))
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
 END DO
 CLOSE(15)
 ! calculation of the stellar radius and Rinf
 min_radius = 1.D99
 DO I = 1, n_modelgrid
  act_radius = model_grid(I)%rwind
  IF(act_radius < min_radius) min_radius = act_radius
 END DO
 R_star = min_radius
 ! R_inf
 max_radius = 1.D0
 DO I = 1, n_modelgrid
  act_radius = model_grid(I)%rwind
  IF(act_radius > max_radius) THEN
   max_radius = act_radius
   max_radius_index = I
  END IF
 END DO
 max_z = 1.D0
 DO I = 1, n_modelgrid
  act_z = model_grid(I)%zwind
  IF(act_z > max_z) THEN
   max_z = act_z
  END IF
 END DO
 write(99,*) 'R_star = ', R_star
 V_inf = model_grid(max_radius_index)%vel
 R_star = MINVAL(model_grid(:)%rwind)
 R_inf = MAXVAL(model_grid(:)%rwind)
 Z_inf = R_inf
 write(*,*) 'read_2D_basic: R_star = ', R_star, ' R_inf = ', R_inf, 'R_inf/R_star = ', R_inf/R_star
 write(99,*) 'computed R_star = ', R_star, ' R_inf = ', R_inf

 xmax = R_inf
 ymax = R_inf
 zmax = R_inf
 xmin = -xmax
 ymin = -ymax
 zmin = -zmax

 ! setting up vacuum and outward model cells
 model_grid(n_modelgrid + 1)%assoc_cells = 0
 model_grid(n_modelgrid + 1)%rwind = 0.D0
 model_grid(n_modelgrid + 1)%zwind = 0.D0
 model_grid(n_modelgrid + 1)%vel   = 0.D0
 model_grid(n_modelgrid + 1)%velang   = 0.D0
 model_grid(n_modelgrid + 1)%rho   = 0.D0
 ! vacuum grids
 model_grid(n_modelgrid + add_mg)%assoc_cells = 0
 model_grid(n_modelgrid + add_mg)%rwind = 0.D0
 model_grid(n_modelgrid + add_mg)%zwind = 0.D0
 model_grid(n_modelgrid + add_mg)%vel   = 0.D0
 model_grid(n_modelgrid + add_mg)%velang   = 0.D0
 model_grid(n_modelgrid + add_mg)%rho   = 0.D0
 model_grid(n_modelgrid + add_mg)%T = 0.D0
 model_grid(n_modelgrid + add_mg)%J= 0.D0
 ALLOCATE (model_grid(n_modelgrid + add_mg)%grid_comp(n_elements))
  DO J = 1, n_elements
   numbions = elements(J)%nions
   ALLOCATE (model_grid(n_modelgrid + add_mg)%grid_comp(J)%grid_ion(numbions))
   atom_number = elements(J)%atom_number
   model_grid(n_modelgrid + 2)%grid_comp(J)%abund = 0.D0
   ! Calculate total number density for included species
   ! tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
   ! model_grid(I)%grid_comp(J)%numb_den = tot_nd
  END DO

 ! OPEN(77, FILE='input_2D_model.dat')
 !  DO cur_mgi = 1, n_modelgrid
 !   write(77,*) model_grid(cur_mgi)%rwind, model_grid(cur_mgi)%angle, model_grid(cur_mgi)%rho, model_grid(cur_mgi)%T, &
 !    model_grid(cur_mgi)%vel, model_grid(cur_mgi)%velang
 !  END DO
 ! CLOSE(77)

 END SUBROUTINE read_2D_basic
