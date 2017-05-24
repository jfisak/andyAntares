! this subroutine calculates collisional rates for the given energy
! level
SUBROUTINE collisional_rates(approx, pack_index, level, nlns, linetransitions, nluns, &
lineuptransitions, population, Ldown, Zdown, Lup, Zup, Lcoll, Zcoll)
USE types
IMPLICIT NONE
! input variables
! used approximation for the collisional term calculation
INTEGER                                 :: approx
! line -- number of line in the linelist field
INTEGER                                 :: pack_index, level, nlns, nluns
INTEGER, DIMENSION(nlns)                :: linetransitions, lineuptransitions
DOUBLE PRECISION                        :: population
! constans
DOUBLE PRECISION, PARAMETER             :: c0 = 5.465D-11
DOUBLE PRECISION, PARAMETER             :: IH = 13.6 * e_v
DOUBLE PRECISION, PARAMETER             :: coll_const = 14.5
! indexes
INTEGER                                 :: act_line
! atomic data
INTEGER                                 :: element_index, ion_index
! value of collision rate
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: current_mgi
! physical parameters of the cell
DOUBLE PRECISION                        :: el_temperature, electron_density, &
                                           temperature
INTEGER                                 :: get_package_model_index
! loop variable
INTEGER                                 :: I
! a value of gamma function
DOUBLE PRECISION                        :: gf
! parameters of transitions
DOUBLE PRECISION                        :: freq, osc_str
DOUBLE PRECISION                        :: exci_energy_l, exci_energy_u
DOUBLE PRECISION                        :: x
! output variables
! total rate
DOUBLE PRECISION                        :: Zdown, Zup, Zcoll
! the rates for the given transitions
DOUBLE PRECISION, DIMENSION(nlns)       :: Ldown, Lcoll
DOUBLE PRECISION, DIMENSION(nluns)      :: Lup
DOUBLE PRECISION                        :: stat_weight


SELECT CASE(approx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! van Regemorter approximation
CASE(1)
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 electron_density = model_grid(current_mgi)%e_dens
 ! --
 temperature = model_grid(current_mgi)%T
 el_temperature = temperature
 ! number of lines we are interested in
 nlns = SIZE(linetransitions)
 element_index = linelist(linetransitions(1))%indexe
 ion_index = linelist(linetransitions(1))%indexi
 Zup = 0.D0
 Zdown = 0.D0
 Zcoll = 0.D0
 DO I = 1, nlns
  ! important physical quantities
  act_line = linetransitions(I)
  osc_str = linelist(act_line)%f_ul
  stat_weight = elements(element_index)%ions(ion_index)%levels(level)%stat_waight
  exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
  exci_energy_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
  ! frequency of transition
  freq = linelist(act_line)%freq
  x = (h * freq) / (BOLK * temperature)
  ! gamma function
  CALL gamma_function(x, act_line, gf)
  ! value of the collision coefficient
  actVal = population * electron_density * c0 * (temperature)**(1.0/2.0) * &
        coll_const * (IH / (h * freq)) * osc_str * &
        ((h * freq) / (BOLK * el_temperature)) * &
        exp(-(h * freq) / (BOLK * el_temperature)) * gf
!  actVal = actVal * (linelist(act_line)%upper - linelist(act_line)%lower)
 ! internal downward jump
  Ldown(I) = actVal * exci_energy_l
  Zdown = Zdown + Ldown(I) * stat_weight
 ! collisional deexcitation
  Lcoll(I) = actVal * (exci_energy_u - exci_energy_l)
  Zcoll = Zcoll + Lcoll(I) * stat_weight
!        (linelist(act_line)%upper - linelist(act_line)%lower)
! print*, 'cool_excit: pop = ', pop, ' electron_density = ', electron_density, &
!        ' temperature = ', temperature, ' gf = ', gf, ' osc_str = ', osc_str, &
!        ' x = ', x
!  print*, 'collisional_rates: actVal = ', actVal
 END DO
! print*, 'collisional_rates: Ztot/pop = ', Ztot
! Ztot = pop * Ztot
! print*, 'collisional_rates: Ztot = ', Ztot
 ! upward jumps
 DO I = 1, nluns
  act_line = lineuptransitions(I)
  ! important physical quantities
  osc_str = linelist(linetransitions(act_line))%f_ul
  stat_weight = elements(element_index)%ions(ion_index)%levels(level)%stat_waight
  exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
  ! frequency of transition
  freq = linelist(linetransitions(act_line))%freq
  x = (h * freq) / (BOLK * temperature)
  ! gamma function
  CALL gamma_function(x, act_line, gf)
  ! value of the collision coefficient
  actVal = population * electron_density * c0 * (temperature)**(1.0/2.0) * &
        coll_const * (IH / (h * freq)) * osc_str * &
        ((h * freq) / (BOLK * el_temperature)) * &
        exp(-(h * freq) / (BOLK * el_temperature)) * gf
  Lup(I) = actVal * exci_energy_l
  Zup = Zup + Lup(I) * stat_weight
 END DO
CASE DEFAULT
 STOP 'collisional_rates: this approximation is not known'
END SELECT


END SUBROUTINE collisional_rates
