SUBROUTINE save_rates()

USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: cur_line
INTEGER                                 :: indexe, indexi
INTEGER                                 :: cur_mgi
INTEGER                                 :: lower_level, upper_level
DOUBLE PRECISION                        :: stat_weight_l, stat_weight_u
DOUBLE PRECISION                        :: ROverV, corrFactor, cur_kappa, cur_radius, cur_rho
DOUBLE PRECISION                        :: constanta
DOUBLE PRECISION                        :: cur_wavelength
DOUBLE PRECISION                        :: exci_energy_u, exci_energy_l
DOUBLE PRECISION                        :: low_pop, upp_pop
DOUBLE PRECISION                        :: roverw

CHARACTER(LEN=120)                      :: file_linetrans

DOUBLE PRECISION                        :: taulu, betalu, Aul, Jlu
DOUBLE PRECISION                        :: actVal
DOUBLE PRECISION                        :: i_downrad, i_int_down, i_int_up
DOUBLE PRECISION                        :: Blu, Bul, fr_line

DOUBLE PRECISION                        :: flux_function

constanta = (const_pi * e_charge**2)/( const_me_g * light_speed)

file_linetrans = trim(outputfolder)//'/r_linetrans.dat'

OPEN(101,FILE=file_linetrans)
write(101,*) '* radius, wavelength, tau_line, i_radeexc, i_int_down, i_int_up'
write(101,*) n_modelgrid
DO cur_mgi = 1, n_modelgrid

  cur_radius = model_grid(cur_mgi)%rwind
  cur_rho = model_grid(cur_mgi)%rho
  ! r opacities
  
  ! line optical depths
 DO cur_line = 1, ntransitions
  fr_line = linelist(cur_line)%freq
  cur_wavelength = 1e8 * light_speed / fr_line
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

  cur_kappa = constanta * linelist(cur_line)%f_lu * low_pop * corrFactor / cur_rho


  ! i opacities
  taulu = light_speed / fr_line * constanta * &
   linelist(cur_line)%f_lu * upp_pop * ROverV * corrFactor
  betalu = 1.D0 / taulu * (1.D0 - exp(- taulu))
  Aul = 8.D0 * fr_line**2 * const_pi**2 * e_charge**2/ (const_me_g * light_speed**3) *&
   stat_weight_l / stat_weight_u * linelist(cur_line)%f_lu
  Jlu = flux_function(0, fr_line, model_grid(cur_mgi)%T, model_grid(cur_mgi)%rwind)

  actVal = Aul * betalu * upp_pop
  
  i_downrad = actVal * exci_energy_l
  i_int_down = actVal * (exci_energy_u - exci_energy_l)

  ! internal upward jump rate

  
  Blu = 4 * const_pi**2 * e_charge**2 / (const_me_g * light_speed * const_h * fr_line) * linelist(cur_line)%f_lu
  Bul = DBLE(stat_weight_l) / DBLE(stat_weight_u) * Blu

  actVal = (Blu * low_pop - Bul * upp_pop) * betalu * Jlu

  i_int_up = actVal * exci_energy_l

  write(101,*) cur_radius, cur_wavelength, cur_kappa, i_downrad, i_int_down, i_int_up
 
 END DO
END DO
CLOSE(101)













END SUBROUTINE save_rates
