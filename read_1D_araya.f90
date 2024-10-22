SUBROUTINE read_1D_araya()

! file format
! Teff
! R_star
! R/R_* -- V/km/s -- -- rho/g/cm^3 -- -- -- -- -- -- -- -- -- -- --

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: n_mg_points

INTEGER                                 :: cur_line, reading_grid, J
INTEGER                                 :: numbions, atom_number

DOUBLE PRECISION                        :: effective_temperature, stellar_radius
DOUBLE PRECISION                        :: radius, density, djunk
DOUBLE PRECISION                        :: velocity

CHARACTER(LEN=100)                      :: junk

add_mg = 1

n_mg_points = 0
OPEN(38, FILE=inputmodelFile)

 READ(38,*) junk
 READ(38,*) junk
 ! number of model points
 DO 
  READ(38,*,iostat=reading_grid) junk
  IF(reading_grid /= 0) EXIT
  n_mg_points = n_mg_points + 1
 END DO

 ALLOCATE(model_grid(n_mg_points + add_mg))
 n_modelgrid = n_mg_points
 ! write(*,*) 'read_1D_araya: n_mg_points = ', n_mg_points

 REWIND(38)

 READ(38,*) effective_temperature
 READ(38,*) stellar_radius

 T_eff = effective_temperature
 R_star = stellar_radius
 ! write(*,*) 'read_1D_model: R_star = ', R_star/const_Rsun

 DO cur_line = 1, n_mg_points
  READ(38,*) radius, djunk, velocity, djunk, djunk, density, junk
  model_grid(cur_line)%rwind = radius * R_star
  write(*,*) 'read_1D_araya: r = ', radius * R_star
  model_grid(cur_line)%vel = velocity * 1.D5 ! [velocity] = km/s
  model_grid(cur_line)%rho = density
  model_grid(cur_line)%T = T_eff
  ! write(*,*) 'read_1D_model: r = ', radius * R_star / const_Rsun
  ALLOCATE (model_grid(cur_line)%grid_comp(n_elements))
  DO J = 1, n_elements      
   numbions = elements(J)%nions
   ! write(*,*) 'read_1D_model: numbions = ', numbions
   ALLOCATE (model_grid(cur_line)%grid_comp(J)%grid_ion(numbions))
   atom_number = elements(J)%atom_number
   !model_grid(cur_line)%grid_comp(J)%abund = massfrac(atom_number)        
   model_grid(cur_line)%grid_comp(J)%abund = elements(J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(cur_line)%grid_comp(J)%abund / elements(J)%atom_mass 
   !model_grid(cur_line)%grid_comp(J)%numb_den = tot_nd
  END DO
 END DO
 R_inf = model_grid(n_mg_points)%rwind
 V_inf = model_grid(n_mg_points)%vel * 100.0
 model_grid(n_modelgrid+add_mg)%rwind = 0.D0
 model_grid(n_modelgrid+add_mg)%vel   = 0.D0
 model_grid(n_modelgrid+add_mg)%rho   = 0.D0     
 write(*,*) 'read_1D_araya: R_inf = ', R_inf/R_star, ' V_inf = ', V_inf

 ! saving n_mg_points into the global variable n_modelgrid
 n_modelgrid = n_mg_points







CLOSE(38)

END SUBROUTINE read_1D_araya
