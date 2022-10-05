SUBROUTINE save_rates()

USE types
IMPLICIT NONE

INTEGER                                 :: cur_line
INTEGER                                 :: indexe, indexi
INTEGER                                 :: cur_mgi
INTEGER                                 :: lower_level, upper_level, stat_weight_u
INTEGER                                 :: stat_weight_l
DOUBLE PRECISION                        :: ROverV, corrFactor, cur_tau, cur_radius
DOUBLE PRECISION                        :: constanta
DOUBLE PRECISION                        :: cur_wavelength
DOUBLE PRECISION                        :: exci_energy_u, exci_energy_l
DOUBLE PRECISION                        :: low_pop, upp_pop
DOUBLE PRECISION                        :: roverw

CHARACTER(LEN=120)                      :: file_r_linetrans

constanta = (pi * e_charge**2)/( me_g * light_speed)

file_r_linetrans = trim(outputfolder)//'/r_linetrans.dat'

OPEN(101,FILE=file_r_linetrans)
DO cur_mgi = 1, n_modelgrid

  cur_radius = model_grid(cur_mgi)%rwind
  ! r opacities
  
  ! line optical depths
 DO cur_line = 1, ntransitions
  cur_wavelength = 1e8 * light_speed / linelist(cur_line)%freq
  indexe = linelist(cur_line)%indexe
  indexi = linelist(cur_line)%indexi
  lower_level = linelist(cur_line)%lower
  upper_level = linelist(cur_line)%upper
  stat_weight_u = &
   elements(indexe)%ions(indexi)%levels(upper_level)%stat_waight
  stat_weight_l = &
   elements(indexe)%ions(indexi)%levels(lower_level)%stat_waight
  exci_energy_u = &
   elements(indexe)%ions(indexi)%levels(upper_level)%exci_energy
  exci_energy_l = &
   elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy
  ROverV = roverw()
  CALL populations(indexe, indexi, lower_level, cur_mgi, low_pop)
  CALL populations(indexe, indexi, upper_level, cur_mgi, upp_pop)

  corrFactor = 1.D0 - (stat_weight_l * upp_pop) / (stat_weight_u * low_pop)

  cur_tau = light_speed / linelist(cur_line)%freq * constanta * &
   linelist(cur_line)%f_lu * low_pop * corrFactor * ROverV

  write(101,*) cur_radius, cur_wavelength, cur_tau
 
 END DO
END DO
CLOSE(101)













END SUBROUTINE save_rates
