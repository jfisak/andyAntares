! collisional ionization rates
! 
! INPUT: approximation(INT): obvious
!        pack_index(INT): packet index
! OUTPUT: Zion(DBLE): total cooling ionization index
!         actikrates(krates): rates for transitions
! 
SUBROUTINE cool_ionization(approximation, pack_index, Zion, actikrates)
USE types
USE constants
USE rates_k
IMPLICIT NONE

INTEGER                         :: approximation, pack_index
DOUBLE PRECISION                :: act_pop
INTEGER                         :: indexe, indexi, indexl
! informations about ions and levels
INTEGER                         :: n_ions, n_levels
! grid informations
INTEGER                         :: current_mgi
DOUBLE PRECISION                :: el_dens, temp, factor_x
INTEGER                         :: nfreq
DOUBLE PRECISION                :: gindex
INTEGER                         :: get_package_model_index
INTEGER                         :: ind_I, act_rate
DOUBLE PRECISION, PARAMETER     :: coll_const = 1.55D13
DOUBLE PRECISION, ALLOCATABLE   :: crossfreq(:)
DOUBLE PRECISION                :: eif, freq
! linear interpolation
DOUBLE PRECISION                :: ali, bli, freq1, freq2, func1, func2
INTEGER                         :: actPoint
DOUBLE PRECISION                :: cross_sect
! output variables
DOUBLE PRECISION                :: Zion
TYPE(krates)                    :: actikrates

  
SELECT CASE(approximation)
!_______________________________________________________________
! van regemorter approximation
! cross section from file saved as a table
CASE (1)
 ! number of computed rates
 act_rate = 0
 ! physical informations in the propagation cell
 ! actual model grid index
 current_mgi = get_package_model_index(pack_index)
 ! electron density
 el_dens = model_grid(current_mgi)%e_dens
 ! temperature
 temp = model_grid(current_mgi)%T
 Zion = 0.D0
 DO indexe = 1, n_elements
  n_ions = SIZE(elements(indexe)%ions)
  DO indexi = 1, n_ions - 1
   n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
   DO indexl = 1, n_levels
    !write(*,*) 'cool_ionization: indexe = ', indexe, ' indexi = ', indexi, ' indexl = ', indexl, &
    ! ' act_rate = ', act_rate
    ! number of points in the photcross array
    IF(ALLOCATED(elements(indexe)%ions(indexi)%levels(indexl)%photcros)) THEN
     nfreq = SIZE(elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:))
    ELSE
     nfreq = 0
    END IF
    IF(nfreq /= 0) THEN
     ALLOCATE(crossfreq(nfreq))
     ! number density of a ground state of ion indexi + 1, indexe
     ! frequency
     freq = (MINVAL(elements(indexe)%ions(indexi + 1)%levels(:)%exci_energy) - &
              elements(indexe)%ions(indexi)%levels(indexl)%exci_energy) / const_h
     !write(*,*) 'cool_ionization: indexe = ', indexe, ' indexi = ', indexi, ' indexl = ', indexl, ' freq = ', freq
     ! argument of E_1(x)
     factor_x = (const_h * freq) / (const_kB * temp)
     ! exponential integral function (calculation of eif)
     CALL exp_int_func(1, factor_x, eif)
     ! write(*,*) 'cool_ionization: calling populations...'
     CALL populations(indexe, indexi, indexl, current_mgi, act_pop)
     ! photoionization cross section
     !CALL bound_free_rates(pack_index, indexl, rad_rate)
     ! photoionization cross section for the given frequency 
     crossfreq(:) = elements(indexe)%ions(indexi)%levels(indexl)%photcros(1,:)
     actPoint = 0
     DO ind_I = 1, nfreq
      IF(freq < crossfreq(ind_I)) THEN
       !write(*,*) 'cool_ionization: freq = ', freq, ' crossfreq(', I, ') = ', crossfreq(I)
       actPoint = ind_I
       EXIT
      END IF
     END DO
     ! if the frequency is out of range of the frequency interval
     ! the total rate will be equal to zero
     IF(actPoint == 0 .OR. actPoint == 1) THEN
      act_rate = act_rate + 1
      actikrates%Lcool_ion(act_rate) = 0.D0
     ELSE
      ! now we have to do a linear interpolation between the points actPoint - 1 and actPoint
      freq1 = crossfreq(actPoint - 1)
      freq2 = crossfreq(actPoint)
      func1 = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2, actPoint - 1)
      func2 = elements(indexe)%ions(indexi)%levels(indexl)%photcros(2, actPoint)
      ali = (func1 - func2) / (freq1 - freq2)
      bli = (func2 * freq1 - func1 * freq2) / (freq1 - freq2)
      cross_sect = ali * freq + bli
      !print*, 'photoionization: ali = ', ali, ' bli = ', bli, ' cross_sect = ', cross_sect
      ! gindex
      IF(indexi == 1) THEN
       gindex = 1.D-1
      ELSE IF(indexi == 2) THEN
       gindex = 2.D-1
      ELSE IF(indexi > 2) THEN
       gindex = 3.D-1
      END IF
      act_rate = act_rate + 1
      actikrates%Lcool_ion(act_rate) = act_pop * el_dens * coll_const / temp**(1.0/2.0) * gindex * &
       cross_sect * exp(-factor_x) / factor_x * (const_h * freq) 
      !write(*,*) 'cool_ionization: act_pop = ', act_pop, ' el_dens = ', el_dens, ' temp = ', temp,&
      ! ' cross_sect = ', cross_sect
      !write(*,*) 'cool_ionization: actikrates%Lcool_ion(', act_rate, ') = ', actikrates%Lcool_ion(act_rate)
      Zion = Zion + actikrates%Lcool_ion(act_rate)
      !print*, 'collion_rates: Zion = ', Zion
     END IF ! actPoint == 0 or actPoint == 1
     DEALLOCATE(crossfreq)
    ELSE ! n_freq == 0
     act_rate = act_rate + 1
     actikrates%Lcool_ion(act_rate) = 0.D0
    END IF ! n_freq /= 0
    !write(*,*) 'cool_ionization: actikrates%Lcool_ion(', act_rate, ') = ', actikrates%Lcool_ion(act_rate)
   END DO ! levels
  END DO ! ions
 END DO ! elements
 
CASE DEFAULT
END SELECT


END SUBROUTINE cool_ionization
