SUBROUTINE cool_excit(approx, pack_index, Ztot, Lcoll)
USE types
IMPLICIT NONE
! input variables
INTEGER                                         :: approx, pack_index
! mgi informations
DOUBLE PRECISION                                :: electron_density, temperature, &
                                                   el_temperature, pop
INTEGER                                         :: get_package_model_index
! line information
INTEGER                                         :: line, element_index, ion_index, &
                                                   current_mgi
DOUBLE PRECISION                                :: osc_str, freq 
! radiative rates
DOUBLE PRECISION                                :: actVal, x, gf
! constans
DOUBLE PRECISION, PARAMETER             :: c0 = 5.465D-11
DOUBLE PRECISION, PARAMETER             :: IH = 13.6 * e_v
DOUBLE PRECISION, PARAMETER             :: coll_const = 14.5
! output variables
DOUBLE PRECISION                                :: Ztot
DOUBLE PRECISION, DIMENSION(ntransitions)       :: Lcoll
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
! print*, 'cool_excit: current_mgi = ', current_mgi

DO line = 1, ntransitions
 ! population of the given state
 element_index = linelist(line)%indexe
 ion_index = linelist(line)%indexi
 CALL populations(element_index, ion_index, linelist(line)%upper, current_mgi, pop)
 ! oscilator strength
 osc_str = linelist(line)%f_ul
 ! frequency of transition
 freq = linelist(line)%freq
 x = (h * freq) / (BOLK * temperature)
 ! gamma function
 CALL gamma_function(x, line, gf)
 
 actVal = pop * electron_density * c0 * (temperature)**(1.0/2.0) * &
       coll_const * (IH / (h * freq)) * osc_str * &
       ((h * freq) / (BOLK * el_temperature)) * &
       exp(-(h * freq) / (BOLK * el_temperature)) * gf
! print*, 'cool_excit: pop = ', pop, ' electron_density = ', electron_density, &
!        ' temperature = ', temperature, ' gf = ', gf, ' osc_str = ', osc_str, &
!        ' x = ', x
 Lcoll(line) = actVal
 Ztot = Ztot + actVal
END DO

CASE DEFAULT
 STOP 'cooling rates: this approximation is not known'
END SELECT

END SUBROUTINE
