! this subroutine calculates collisional rates for the given energy
! level
SUBROUTINE collisional_rates(approx, pack_index, act_state, nlns, linetransitions, Ztot, Lcoll)
USE types
IMPLICIT NONE
! input variables
INTEGER                                 :: approx
INTEGER                                 :: pack_index, act_state, nlns
INTEGER, DIMENSION(nlns)                :: linetransitions
! constans
DOUBLE PRECISION, PARAMETER             :: c0 = 5.465D-11
DOUBLE PRECISION, PARAMETER             :: IH = 13.6 * e_v
DOUBLE PRECISION, PARAMETER             :: coll_const = 14.5
! indexes
INTEGER                                 :: line
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
DOUBLE PRECISION                        :: Ztot
! the rates for the given transitions
DOUBLE PRECISION, DIMENSION(nlns) :: Lcoll


SELECT CASE(approx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! van Regemorter approximation
CASE(0)
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 electron_density = model_grid(current_mgi)%e_dens
 ! --
 temperature = model_grid(current_mgi)%T
 el_temperature = temperature
 ! number of lines we are interested in
 nlns = SIZE(linetransitions)

 DO line = 1, nlns
  ! oscilator strength
  osc_str = linelist(linetransitions(line))%f_ul
  ! frequency of transition
  freq = linelist(linetransitions(line))%freq
  x = (h * freq) / (BOLK * temperature)
  ! gamma function
  CALL gamma_function(x, act_state, gf)
  ! value of the collision coefficient
  actVal = electron_density * c0 * (temperature)**(1.0/2.0) * &
        coll_const * (IH / (h * freq)) * osc_str * &
        ((h * freq) / (BOLK * el_temperature)) * &
        exp(-(h * freq) / (BOLK * el_temperature)) * gf
  Lcoll(line) = actVal
  Ztot = Ztot + actVal
 END DO
CASE DEFAULT
 STOP 'collisional_rates: this approximation is not known'
END SELECT


END SUBROUTINE collisional_rates
