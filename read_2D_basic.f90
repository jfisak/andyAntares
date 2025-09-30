! reads a basic 2D model
!
! input file
!
! radius, perpend, dens, velrad, velang, temp
! radius: r * unit_length * const_Rsun
! perpend: ?
! density: rho * unit_density
! velrad: v * unit_velocity
! velang: v * unit_velocity
! temp: T
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE read_2D_basic()

USE types
USE constants


IMPLICIT NONE

INTEGER                                :: ios
INTEGER, PARAMETER                     :: maxrows = 60000000
DOUBLE PRECISION                          :: junk
! loop variables
INTEGER                                :: ind_I, ind_J
! readen physical quantities
DOUBLE PRECISION                       :: radius, perpend, dens, velrad, velang, temp
INTEGER                                :: atom_number, numbions
! calcultaion a stellar radius
DOUBLE PRECISION                       :: act_radius, min_radius, max_radius
INTEGER                                :: max_radius_index

DOUBLE PRECISION                        :: unit_length, unit_velocity, unit_density
DOUBLE PRECISION                        :: min_rad

unit_length = 12.64759321736591 * const_Rsun
unit_velocity = 1.D8
unit_density = 1.41314878888978971775872825028550241D-0006



 add_mg = 3
 write(99,*) 'we will read input input data from the basic 2D model'
 ! firstly we calculate number of rows in the file
 n_modelgrid = 0
 T_eff = 30000
 min_radius = 1.D99
 OPEN(UNIT=15,status='old', FILE=inputmodelFile)
  DO ind_I = 1, maxrows
   READ(15,*,IOSTAT = ios) radius, junk, junk, junk, junk, junk
   if(ios /= 0) EXIT
   radius = radius * unit_length
   IF(radius < min_radius) THEN
    min_radius = radius
   END IF
  END DO
  R_star = min_radius
  min_rad = 3 * R_star
  write(*,*) 'read_2D_basic: R_star = ', R_star, ' min_rad = ', min_rad
  REWIND(15)
  DO ind_I = 1, maxrows
   READ(15,*,IOSTAT = ios) radius, junk, junk, junk, junk, junk
   if(ios /= 0) EXIT
   if(ind_I == maxrows) THEN
    write(99,*) 'maximum number of records exceeded in subroutine read_2d_model'
    write(99,*) 'exiting program now...'
    STOP
   end if
  IF(radius * unit_length > min_rad) THEN
   n_modelgrid = n_modelgrid + 1
  END IF
 END DO
 write(99,*) 'mumber of model grids: ', n_modelgrid
 write(*,*) 'mumber of model grids: ', n_modelgrid
 IF (n_modelgrid .EQ. 0) STOP 'no model grid cells were found...'
 ! n_modelgrid + 1 ... for dummy cells
 ! n_modelgrid + 2 ... for cells with R_star < r < R_inf but too far from some model grid point
 !                     (vacuum cell) 
 ALLOCATE(model_grid(n_modelgrid + add_mg))
 ! define special indeces
 outerspace_index = n_modelgrid + 1
 photosphere_index = n_modelgrid + 2
 vacuum_index = n_modelgrid + 3

 REWIND(15)
 ind_I = 0
 DO 
  READ(15,*, IOSTAT = ios) radius, perpend, dens, velrad, velang, temp
  if(ios /= 0) EXIT
  ! we want a model not starting so deep in photosphere
  IF(radius * unit_length > min_rad) THEN ! if1
   ind_I = ind_I + 1
   model_grid(ind_I)%rwind = radius * unit_length
   model_grid(ind_I)%angle = perpend
   model_grid(ind_I)%vel = velrad * unit_velocity
   model_grid(ind_I)%velang = velang * unit_velocity
   model_grid(ind_I)%rho = dens * unit_density
   model_grid(ind_I)%T = temp
   model_grid(ind_I)%J = 0.D0
   model_grid(ind_I)%assoc_cells = 0
   IF(temp < 1000) STOP 'read_2D_basic: temp < 1000'
   ALLOCATE (model_grid(ind_I)%grid_comp(n_elements))
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
  END IF ! if1: r * n_r > min_r
 END DO
 CLOSE(15)
 ! calculation of the stellar radius and Rinf
 ! R_inf
 max_radius = 1.D0
 DO ind_I = 1, n_modelgrid
  act_radius = model_grid(ind_I)%rwind
  IF(act_radius > max_radius) THEN
   max_radius = act_radius
   max_radius_index = ind_I
  END IF
 END DO
 R_inf = max_radius
!  max_z = 1.D0
!  DO ind_I = 1, n_modelgrid
!   act_z = model_grid(ind_I)%zwind
!   IF(act_z > max_z) THEN
!    max_z = act_z
!   END IF
!  END DO
!  Z_inf = max_z
!  write(99,*) 'R_star = ', R_star
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
 ! photosphere index
 model_grid(photosphere_index)%assoc_cells = 0
 model_grid(photosphere_index)%rwind = 0.D0
 model_grid(photosphere_index)%zwind = 0.D0
 model_grid(photosphere_index)%vel   = 0.D0
 model_grid(photosphere_index)%velang   = 0.D0
 model_grid(photosphere_index)%rho   = 0.D0
 model_grid(outerspace_index)%assoc_cells = 0
 model_grid(outerspace_index)%rwind = 0.D0
 model_grid(outerspace_index)%zwind = 0.D0
 model_grid(outerspace_index)%vel   = 0.D0
 model_grid(outerspace_index)%velang   = 0.D0
 model_grid(outerspace_index)%rho   = 0.D0
 ! vacuum grids
 model_grid(vacuum_index)%assoc_cells = 0
 model_grid(vacuum_index)%rwind  = 0.D0
 model_grid(vacuum_index)%zwind  = 0.D0
 model_grid(vacuum_index)%vel    = 0.D0
 model_grid(vacuum_index)%velang = 0.D0
 model_grid(vacuum_index)%rho    = 0.D0
 model_grid(vacuum_index)%T      = 0.D0
 model_grid(vacuum_index)%J      = 0.D0
 ALLOCATE (model_grid(vacuum_index)%grid_comp(n_elements))
  DO ind_J = 1, n_elements
   numbions = elements(ind_J)%nions
   ALLOCATE (model_grid(n_modelgrid + add_mg)%grid_comp(ind_J)%grid_ion(numbions))
   atom_number = elements(ind_J)%atom_number
   model_grid(vacuum_index)%grid_comp(ind_J)%abund = 0.D0
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
