! calculates rates for the ionic transitions
! 
! INPUT: indexe(INT): element index
!        indexi(INT): ion index
!        leveli(INT): level index
!        current_mgi(INT): modGrid index
!        act_pop(DBLE): current population
! OUTPUT: Zion(DBLE): total ion rate
!         Zintrecom(DBLE): total internal recombination rate
!         Zrecom(DBLE): total recombination rate
!         actirates(DBLE): rates for corresponding transitions
!
SUBROUTINE i_radion(indexe, indexi, leveli, current_mgi, act_pop, Zion, Zintrecom, Zrecom, actirates)
USE types
USE constants
USE rates_i
IMPLICIT NONE

! input variables
INTEGER                                 :: indexe, indexi, leveli
INTEGER                                 :: current_mgi, nrecom
DOUBLE PRECISION                        :: act_pop
! computing fields
INTEGER                                 :: npoints
INTEGER                                 :: ind_I, ind_K
DOUBLE PRECISION                        :: photRate
DOUBLE PRECISION                        :: temp
DOUBLE PRECISION                        :: phot_cross
DOUBLE PRECISION                        :: exci_energy, gr_exci_energy
DOUBLE PRECISION                        :: pop_number
DOUBLE PRECISION                        :: el_dens
DOUBLE PRECISION                        :: actVal
! output variables
DOUBLE PRECISION                        :: Zion, Zrecom, Zintrecom
TYPE(irates)                           :: actirates
! linear interpolation
INTEGER                                 :: act_index, temp_i
DOUBLE PRECISION                        :: ali, bli, func1, func2
DOUBLE PRECISION                        :: temp1, temp2, Tmax, Tmin
INTEGER                                 :: Ntpoints



! write(*,*) 'i_radion: indexe = ', indexe, ' indexi = ', indexi, ' leveli = ', leveli
el_dens = model_grid(current_mgi)%e_dens
temp = model_grid(current_mgi)%t
Ntpoints = SIZE(i_temps)
Tmin = i_temps(1)
Tmax = i_temps(Ntpoints)
act_index = 0
DO ind_I = 1, n_photcrossect
 IF(iints(ind_I)%indexe == indexe .AND. iints(ind_I)%indexi == indexi &
  .AND. iints(ind_I)%indexl == leveli) THEN
  act_index = ind_I
 END IF
END DO
!print*, 'photion_rates: npoints = ', npoints
! there are no data for photoionization cross section available
! the rates are equal to zero
!print*, 'photion_rates: npoints = ', npoints
IF(act_index /= 0) THEN
 ! index of the array
 IF(Ntpoints /= 1) THEN
  IF(temp == Tmax) THEN
   actVal = iints(act_index)%gammaijk(Ntpoints)
  ELSE IF(temp == Tmin) THEN
   actVal = iints(act_index)%gammaijk(1)
  ELSE
   temp_i = FLOOR((temp * DBLE(Ntpoints - 1) - DBLE(Ntpoints) * Tmin + Tmax)/(Tmax - Tmin))
   ! write(*,*) 'i_radion: temp_i = ', temp_i
   ! write(*,*) 'i_radion: temp = ', temp, ' Tmin = ', Tmin, ' Tmax = ', Tmax
   temp1 = i_temps(temp_i)
   temp2 = i_temps(temp_i + 1)
   func1 = iints(act_index)%gammaijk(temp_i)
   func2 = iints(act_index)%gammaijk(temp_i + 1)
   ali = (func2 - func1) / (temp2 - temp1)
   bli = (func1 * temp2 - func2 * temp1) / (temp2 - temp1)
   phot_cross = ali * temp + bli
  END IF
 ELSE IF(Ntpoints == 1) THEN
  phot_cross = iints(act_index)%gammaijk(1)
  ! write(*,*) 'i_radion: act_index = ', act_index, ' phot_cross = ', phot_cross
 END IF
 photRate = act_pop * actVal
 ! write(*,*) 'i_radion: actVal = ', actVal
 ! write(*,*) 'i_radion: indexe = ', indexe, ' indexi - 1 = ', indexi - 1, 'ind_K = ', K
 ! write(*,*) 'i_radion: act_pop = ', act_pop
 IF(photRate < 0.D0) STOP 'i_radion: photRate < 0'
 Zion = photRate * elements(indexe)%ions(indexi)%levels(leveli)%exci_energy
 ! write(*,*) 'i_radion: photRate = ', photRate, ' Zion = ', Zion
ELSE
 Zion = 0.D0
END IF
 ! write(*,*) 'i_radion: Zion = ', Zion
 ! print*, 'photion_rates: phot_cross = ', phot_cross
 !____________________________________________________________________________
 ! recombination
Zintrecom = 0.D0
Zrecom = 0.D0
IF(indexi > 1) THEN
 nrecom = SIZE(actirates%Lma_recrad)
 ! write(*,*) 'i_radion: nrecom = ', nrecom
 DO ind_K = 1, nrecom
  IF(ALLOCATED(elements(indexe)%ions(indexi - 1)%levels(ind_K)%photcros)) THEN
   npoints = SIZE(elements(indexe)%ions(indexi - 1)%levels(ind_K)%photcros(1,:))
  ELSE
   npoints = 0
  END IF
  ! write(*,*) 'i_radion: npoints = ', npoints
  IF(npoints /= 0) THEN  
   ! write(*,*) 'i_radion: calling populations...'
   act_index = elements(indexe)%ions(indexi - 1)%levels(ind_K)%phfreqi
   CALL populations(indexe, indexi, 1, current_mgi, pop_number)
   IF(Ntpoints /= 1) THEN
    IF(temp == Tmax) THEN
     actVal = iints(act_index)%alphaijk(Ntpoints)
    ELSE IF(temp == Tmin) THEN
     actVal = iints(act_index)%alphaijk(1)
    ELSE
     temp_i = FLOOR((temp * DBLE(Ntpoints - 1) - DBLE(Ntpoints) * Tmin + Tmax)/(Tmax - Tmin))
     IF(temp_i <= 0) THEN
      write(*,*) 'i_radion: temp_i = ', temp_i
      write(*,*) ' temp = ', temp, ' Tmin = ', Tmin, ' Tmax = ', Tmax, ' Ntpoints = ', Ntpoints
     END IF
     ! write(*,*) 'i_radion: temp_i = ', temp_i
     ! write(*,*) 'i_radion: temp = ', temp, ' Tmin = ', Tmin, ' Tmax = ', Tmax
     temp1 = i_temps(temp_i)
     temp2 = i_temps(temp_i + 1)
     func1 = iints(act_index)%alphaijk(temp_i)
     func2 = iints(act_index)%alphaijk(temp_i + 1)
     ali = (func2 - func1) / (temp2 - temp1)
     bli = (func1 * temp2 - func2 * temp1) / (temp2 - temp1)
     phot_cross = ali * temp + bli
    END IF
   ELSE IF(Ntpoints == 1) THEN
    phot_cross = iints(act_index)%alphaijk(1)
    ! write(*,*) 'i_radion: act_index = ', act_index, ' phot_cross = ', phot_cross
   END IF
   ! write(*,*) 'i_radion: temp1 ', temp1, ' temp2 = ', temp2, ' func1 = ', func1, &
   !  ' func2 = ', func2
   ! write(*,*) 'i_radion: ali = ', ali, ' bli = ', bli, 'phot_cross = ', phot_cross
   exci_energy = elements(indexe)%ions(indexi - 1)%levels(ind_K)%exci_energy
   gr_exci_energy = MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
   actVal = pop_number * el_dens * phot_cross
   actirates%Lma_int_recrad(ind_K) = actVal * exci_energy
   actirates%Lma_recrad(ind_K) = actVal * (gr_exci_energy - exci_energy)
   !write(*,*) 'i_radion: actVal = ', actVal
   ! write(*,*) 'i_radion: phot_cross = ', phot_cross
   ! write(*,*) 'i_radion: indexe = ', indexe, ' indexi - 1 = ', indexi - 1, 'K = ', K
   ! write(*,*) 'i_radion: exci_energy = ', exci_energy, ' gr_exci_energy = ', gr_exci_energy
   ! write(*,*) 'i_radion: pop_number = ', pop_number, ' el_dens = ', el_dens
   ! write(*,*) 'i_radion: Lint = ', actirates%Lma_int_recrad(ind_K), ' Lrec = ', actirates%Lma_recrad(ind_K)
   Zintrecom = Zintrecom + actirates%Lma_int_recrad(ind_K) 
   Zrecom = Zrecom + actirates%Lma_recrad(ind_K)
   ! print*, 'Zrecom = ', Zrecom
   act_index = act_index + 1
  ELSE ! npoints = 0
   actirates%Lma_int_recrad(ind_K) = 0.D0
   actirates%Lma_recrad(ind_K) = 0.D0
  END IF ! npoints
 END DO
 ! STOP 'i_radion: testing'
END IF ! indexi > 1
! write(*,*)  'photion_rates: Zion = ', Zion, ' Zrecom = ', Zrecom, ' Zintrecom = ', Zintrecom
END SUBROUTINE i_radion
