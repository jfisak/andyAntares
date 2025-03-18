! calculate collisional exciation rates
!
! INPUT: approx(INT): approximation
!        pack_index(INT): index of packet
!        actikrates(krates): contains all rates for all transitions
! OUTPUT: Zexc(INT): total collisional excitational rate
! 
SUBROUTINE cool_excit(approx, pack_index, Zexc, actikrates)
USE types
USE constants
USE rates_k
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
DOUBLE PRECISION                                :: exc_upper, exc_lower
! radiative rates
DOUBLE PRECISION                                :: actVal, factor_x, gf
! constans
DOUBLE PRECISION, PARAMETER                     :: c0 = 5.465D-11
DOUBLE PRECISION, PARAMETER                     :: IH = 13.6 * const_ev
DOUBLE PRECISION, PARAMETER                     :: coll_const = 14.5
! output variables
DOUBLE PRECISION                                :: Zexc
TYPE(krates)                                    :: actikrates

Zexc = 0.D0
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
 ! el_temperature = temperature / 1.D2
 ! print*, 'cool_excit: current_mgi = ', current_mgi
 
 DO line = 1, ntransitions
  ! population of the given state
  element_index = linelist(line)%indexe
  ion_index = linelist(line)%indexi
  ! write(*,*) 'cool_excit: calling populations...'
  CALL populations(element_index, ion_index, linelist(line)%lower, current_mgi, pop)
  ! oscilator strength
  osc_str = linelist(line)%f_lu
  ! frequency of transition
  freq = linelist(line)%freq
  factor_x = (const_h * freq) / (const_kB * temperature)
  exc_upper = elements(element_index)%ions(ion_index)%levels(linelist(line)%upper)%exci_energy
  exc_lower = elements(element_index)%ions(ion_index)%levels(linelist(line)%lower)%exci_energy
  ! gamma function
  CALL gamma_function(factor_x, line, gf)
  
  actVal = electron_density * c0 * (temperature)**(1.0/2.0) * &
        coll_const * (IH / (const_h * freq)) * osc_str * &
        ((const_h * freq) / (const_kB * el_temperature)) * &
        exp(-factor_x) * gf * (exc_upper - exc_lower)
  ! print*, 'cool_excit: pop = ', pop, ' electron_density = ', electron_density, &
  !        ' temperature = ', temperature, ' gf = ', gf, ' osc_str = ', osc_str, &
  !        ' factor_x = ', x
  actikrates%Lcool_excit(line) = pop * actVal
  !write(*,*) 'cool_excit: actikrates%Lcool_excit(', line, ') = ', actikrates%Lcool_excit(line)
  Zexc = Zexc + actVal
  !write(*,*) 'cool_excit: Zexc = ', Zexc
 END DO
CASE DEFAULT
 STOP 'cooling rates: this approximation is not known'
END SELECT

END SUBROUTINE
