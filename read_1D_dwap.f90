! reads 1D model provided by Araya
! file format
! text description of columns
! 0:I 1:r[cm] 2:rho[g/cm3] 3:vr[cm/s] 4:vt[cm/s] 5:Tgas[K] 6:Trad[K]
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE read_1D_dwap()

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: n_mg_points

INTEGER                                 :: cur_line, reading_grid, ind_J
INTEGER                                 :: numbions, atom_number, cur_I

DOUBLE PRECISION                        :: effective_temperature
DOUBLE PRECISION                        :: radius, density, t_gas, t_rad
DOUBLE PRECISION                        :: vel_rad, vel_tan

CHARACTER(LEN=filename_length)                      :: junk

add_mg = 2

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
 write(*,*) 'read_1D_dwap: n_mg_points = ', n_mg_points

 REWIND(38)

 READ(38,*) junk
 READ(38,*) effective_temperature
 
 T_eff = effective_temperature

 DO cur_line = 1, n_mg_points
  READ(38,*) cur_I, radius, density, vel_rad, vel_tan, t_gas, t_rad
  model_grid(cur_line)%rwind = radius
  model_grid(cur_line)%vel = vel_rad ! [velocity] = km/s
  model_grid(cur_line)%velang = vel_tan ! [velocity] = km/s
  model_grid(cur_line)%rho = density
  model_grid(cur_line)%T = t_gas
  ! write(*,*) 'read_1D_model: r = ', radius * R_star / const_Rsun
  ALLOCATE (model_grid(cur_line)%grid_comp(n_elements))
  DO ind_J = 1, n_elements      
   numbions = elements(ind_J)%nions
   ! write(*,*) 'read_1D_model: numbions = ', numbions
   ALLOCATE (model_grid(cur_line)%grid_comp(ind_J)%grid_ion(numbions))
   atom_number = elements(ind_J)%atom_number
   !model_grid(cur_line)%grid_comp(ind_J)%abund = massfrac(atom_number)        
   model_grid(cur_line)%grid_comp(ind_J)%abund = elements(ind_J)%abundance
   !Calculate total number density for included species
   !tot_nd = model_grid(cur_line)%grid_comp(ind_J)%abund / elements(ind_J)%atom_mass 
   !model_grid(cur_line)%grid_comp(ind_J)%numb_den = tot_nd
  END DO
 END DO
 R_star = MINVAL(model_grid(1:n_mg_points)%rwind)
 R_inf = MAXVAL(model_grid(1:n_mg_points)%rwind)
 V_inf = model_grid(n_mg_points)%vel
 outerspace_index = n_modelgrid + 1
 photosphere_index = n_modelgrid + 2
 model_grid(outerspace_index)%rwind = 0.D0
 model_grid(outerspace_index)%vel   = 0.D0
 model_grid(outerspace_index)%rho   = 0.D0     
 model_grid(photosphere_index)%rwind = 0.D0
 model_grid(photosphere_index)%vel   = 0.D0
 model_grid(photosphere_index)%rho   = 0.D0     
 ! write(*,*) 'read_1D_dwap: R_inf = ', R_inf/R_star, ' V_inf = ', V_inf

CLOSE(38)

END SUBROUTINE read_1D_dwap
