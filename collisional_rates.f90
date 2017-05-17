! this subroutine calculates collisional rates for the given energy
! level
SUBROUTINE collisional_rates(approx, lower,pack_index, line, nlns, linetransitions, Ztot, Lcoll)
USE types
IMPLICIT NONE
! input variables
INTEGER                                 :: approx
INTEGER                                 :: pack_index, line, nlns
INTEGER, DIMENSION(nlns)                :: linetransitions
DOUBLE PRECISION                        :: pop
! = 1 for upward = 0 for downward transitions
LOGICAL                                 :: lower
! constans
DOUBLE PRECISION, PARAMETER             :: c0 = 5.465D-11
!DOUBLE PRECISION, PARAMETER             :: c0 = 5.46510**(-11)
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
! a value of gamma function
DOUBLE PRECISION                        :: gf
! parameters of transitions
DOUBLE PRECISION                        :: freq, osc_str
DOUBLE PRECISION                        :: x
! output variables
! total rate
DOUBLE PRECISION                        :: Ztot, stat_weight
! the rates for the given transitions
INTEGER                                 :: level_index
DOUBLE PRECISION, DIMENSION(nlns) :: Lcoll


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
 element_index = linelist(line)%indexe
 ion_index = linelist(line)%indexi
 Ztot = 0.D0
  IF (lower .EQV. .TRUE.) THEN
   level_index = linelist(line)%lower
  ELSE IF (lower .EQV. .FALSE.) THEN
   level_index = linelist(line)%upper
  END IF
 DO act_line = 1, nlns
  CALL populations(element_index, ion_index, level_index, current_mgi, pop)
  ! oscilator strength
  osc_str = linelist(linetransitions(act_line))%f_ul
  stat_weight = elements(element_index)%ions(ion_index)%levels(level_index)%stat_waight
  ! frequency of transition
  freq = linelist(linetransitions(act_line))%freq
  x = (h * freq) / (BOLK * temperature)
  ! gamma function
  CALL gamma_function(x, act_line, gf)
  ! value of the collision coefficient
  actVal = electron_density * c0 * (temperature)**(1.0/2.0) * &
        coll_const * (IH / (h * freq)) * osc_str * &
        ((h * freq) / (BOLK * el_temperature)) * &
        exp(-(h * freq) / (BOLK * el_temperature)) * gf
!  actVal = actVal * (linelist(act_line)%upper - linelist(act_line)%lower)
  Lcoll(act_line) = actVal
  Ztot = Ztot + actVal * stat_weight * &
        (linelist(act_line)%upper - linelist(act_line)%lower)
! print*, 'cool_excit: pop = ', pop, ' electron_density = ', electron_density, &
!        ' temperature = ', temperature, ' gf = ', gf, ' osc_str = ', osc_str, &
!        ' x = ', x
!  print*, 'collisional_rates: actVal = ', actVal
 END DO
 Ztot = pop * Ztot
! print*, 'collisional_rates: Ztot = ', Ztot
CASE DEFAULT
 STOP 'collisional_rates: this approximation is not known'
END SELECT


END SUBROUTINE collisional_rates
