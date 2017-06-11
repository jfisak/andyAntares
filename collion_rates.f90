SUBROUTINE collion_rates(approximation, indexe, indexi, pack_index, act_level, act_pop, Zion)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: approximation
INTEGER                         :: pack_index, act_level
DOUBLE PRECISION                :: act_pop
INTEGER                         :: indexe, indexi
! grid informations
INTEGER                         :: current_mgi
DOUBLE PRECISION                :: el_dens, temp, gl_pop_ip1e, x
INTEGER                         :: gindex, nfreq
INTEGER                         :: get_package_model_index
INTEGER                         :: I
DOUBLE PRECISION, PARAMETER     :: coll_const = 1.55D13
DOUBLE PRECISION, ALLOCATABLE   :: crossfreq(:)
DOUBLE PRECISION                :: eif, freq
! linear interpolation
DOUBLE PRECISION                :: ali, bli, freq1, freq2, func1, func2
INTEGER                         :: actPoint
DOUBLE PRECISION                :: cross_sect
! output variables
DOUBLE PRECISION                :: Zion


SELECT CASE(approximation)

!_______________________________________________________________
! hydrogenic approximation
! cross section is proportional to (f_ijk/f)^3
CASE (1)
 ! actual model grid index
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 el_dens = model_grid(current_mgi)%e_dens
 ! temperature
 temp = model_grid(current_mgi)%T
 ! number density of a ground state of ion indexi + 1, indexe
 gl_pop_ip1e = model_grid(current_mgi)%grid_comp(indexe)%grid_ion(indexi + 1)%gl_pop
 ! frequency
  freq = (elements(indexe)%ions(indexi + 1)%levels(1)%exci_energy - &
          elements(indexe)%ions(indexi)%levels(act_level)%exci_energy) / h
 ! argument of E_1(x)
 x = (h * freq) / (BOLK * temp)
 ! exponential integral function (calculation of eif)
 CALL exp_int_func(1, x, eif)
 ! photoionization cross section
 !CALL bound_free_rates(pack_index, act_level, rad_rate)
 ! photoionization cross section for the given frequency freq
 nfreq = SIZE(elements(indexe)%ions(indexi)%levels(act_level)%photcros(1,:))
 ALLOCATE(crossfreq(nfreq))
 crossfreq(:) = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1,:)
 actPoint = 0
 DO I = 1, nfreq
  IF(freq > crossfreq(I)) THEN
   actPoint = I
   EXIT
  END IF
 END DO
 ! if the frequency is out of range of the frequency interval
 ! the total rate will be equal to zero
 IF(actPoint == 0 .OR. actPoint == 1) THEN
  Zion = 0
  RETURN
 END IF
 ! now we have to do a linear interpolation between the points actPoint - 1 and actPoint
 freq1 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1, I - 1)
 freq2 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(1, I)
 func1 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(2, I - 1)
 func2 = elements(indexe)%ions(indexi)%levels(act_level)%photcros(2, I)
 ali = (func1 - func2) / (freq1 - freq2)
 bli = (func2 * freq1 - func1 * freq2) / (freq1 - freq2)
 cross_sect = ali * freq + bli
 ! gindex
 IF(indexi == 1) THEN
  gindex = 0
 ELSE IF(indexi == 2) THEN
  gindex = 1
 ELSE IF(indexi > 2) THEN
  gindex = 2
 END IF

 Zion = act_pop * el_dens *coll_const / temp**(1.0/2.0) * gindex * cross_sect * exp(-x) / x 
 !print*, 'collion_rates: Zion = ', Zion


CASE DEFAULT
 STOP 'bound_free_rates: non valid approximation was chosen'
END SELECT


END SUBROUTINE collion_rates
