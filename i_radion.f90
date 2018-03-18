SUBROUTINE i_radion(approx, indexe, indexi, leveli, current_mgi, act_pop, Zion, Zintrecom, Zrecom, actirates)
USE types
USE rates_i
IMPLICIT NONE

! input variables
INTEGER                                 :: approx, indexe, indexi, leveli
INTEGER                                 :: current_mgi, nrecom
DOUBLE PRECISION                        :: act_pop
! computing fields
INTEGER                                 :: npoints
INTEGER                                 :: I, J, K
DOUBLE PRECISION, ALLOCATABLE           :: freq(:), cross(:), func(:)
DOUBLE PRECISION                        :: gijk, photRate
DOUBLE PRECISION                        :: temp
DOUBLE PRECISION                        :: up_pop
DOUBLE PRECISION                        :: flux
DOUBLE PRECISION                        :: summ
DOUBLE PRECISION                        :: phot_cross
DOUBLE PRECISION                        :: freqt
DOUBLE PRECISION                        :: exci_energy, gr_exci_energy
DOUBLE PRECISION                        :: pop_number
DOUBLE PRECISION                        :: x
DOUBLE PRECISION                        :: flux_function
DOUBLE PRECISION                        :: stat_weight
DOUBLE PRECISION                        :: el_dens
DOUBLE PRECISION                        :: sfactor
DOUBLE PRECISION                        :: actVal
! output variables
DOUBLE PRECISION                        :: Zion, Zrecom, Zintrecom
TYPE(irates)                           :: actirates
INTEGER                                 :: OMP_GET_THREAD_NUM, my_rank
! linear interpolation
INTEGER                                 :: act_index, temp_i
DOUBLE PRECISION                        :: ali, bli, func1, func2
DOUBLE PRECISION                        :: temp1, temp2, Tmax, Tmin

my_rank = OMP_GET_THREAD_NUM()


! write(*,*) 'i_radion: indexe = ', indexe, ' indexi = ', indexi, ' leveli = ', leveli
el_dens = model_grid(current_mgi)%e_dens
temp = model_grid(current_mgi)%t
Tmin = i_temps(1)
Tmax = i_temps(SIZE(i_temps))
act_index = 0
DO I = 1, n_photcrossect
 IF(iints(I)%indexe == indexe .AND. iints(I)%indexi == indexi &
  .AND. iints(I)%indexl == leveli) THEN
  act_index = I
 END IF
END DO
!print*, 'photion_rates: npoints = ', npoints
! there are no data for photoionization cross section available
! the rates are equal to zero
!print*, 'photion_rates: npoints = ', npoints
IF(act_index /= 0) THEN
 ! index of the array
 temp_i = FLOOR((temp - Tmin) / (Tmax - Tmin) * DBLE(SIZE(i_temps)))
 temp1 = i_temps(temp_i)
 temp2 = i_temps(temp_i + 1)
 func1 = iints(act_index)%gammaijk(temp_i)
 func2 = iints(act_index)%gammaijk(temp_i + 1)
 ali = (func1 - func2) / (temp1 - temp2)
 bli = (func2 * temp1 - func1 * temp2) / (temp1 - temp2)
 actVal = ali * temp + bli
 photRate = act_pop * actVal
 write(*,*) 'i_radion: actVal = ', actVal
 write(*,*) 'i_radion: indexe = ', indexe, ' indexi - 1 = ', indexi - 1, 'K = ', K
 write(*,*) 'i_radion: act_pop = ', act_pop
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
 DO I = 1, n_photcrossect
  IF(iints(I)%indexe == indexe .AND. iints(I)%indexi == indexi - 1) THEN
   act_index = I
   ! write(*,*) 'i_radion: found phcs: ', I
   EXIT
  END IF
 END DO
 nrecom = SIZE(actirates%Lma_recrad)
 ! write(*,*) 'i_radion: nrecom = ', nrecom
 DO K = 1, nrecom
  IF(ALLOCATED(elements(indexe)%ions(indexi - 1)%levels(K)%photcros)) THEN
   npoints = SIZE(elements(indexe)%ions(indexi - 1)%levels(K)%photcros(1,:))
  ELSE
   npoints = 0
  END IF
  ! write(*,*) 'i_radion: npoints = ', npoints
  IF(npoints /= 0) THEN  
   ! write(*,*) 'i_radion: calling populations...'
   CALL populations(indexe, indexi, 1, current_mgi, pop_number)
   temp_i = FLOOR((temp - Tmin) / (Tmax - Tmin) * DBLE(SIZE(i_temps)))
   ! write(*,*) 'i_radion: temp_i = ', temp_i
   temp1 = i_temps(temp_i)
   temp2 = i_temps(temp_i + 1)
   func1 = iints(act_index)%alphaijk(temp_i)
   func2 = iints(act_index)%alphaijk(temp_i + 1)
   ali = (func2 - func1) / (temp2 - temp1)
   bli = (func1 * temp2 - func2 * temp1) / (temp2 - temp1)
   phot_cross = ali * temp + bli
   ! write(*,*) 'i_radion: temp1 ', temp1, ' temp2 = ', temp2, ' func1 = ', func1, &
   !  ' func2 = ', func2
   ! write(*,*) 'i_radion: ali = ', ali, ' bli = ', bli, 'phot_cross = ', phot_cross
   exci_energy = elements(indexe)%ions(indexi - 1)%levels(K)%exci_energy
   gr_exci_energy = MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
   ! actVal = pop_number * el_dens * phot_cross
   actirates%Lma_int_recrad(K) = phot_cross * pop_number * exci_energy
   actirates%Lma_recrad(K) = phot_cross * pop_number * (gr_exci_energy - exci_energy)
   ! write(*,*) 'i_radion: actVal = ', actVal
   ! write(*,*) 'i_radion: indexe = ', indexe, ' indexi - 1 = ', indexi - 1, 'K = ', K
   ! write(*,*) 'i_radion: exci_energy = ', exci_energy, ' gr_exci_energy = ', gr_exci_energy
   ! write(*,*) 'i_radion: pop_number = ', pop_number, ' el_dens = ', el_dens
   ! write(*,*) 'i_radion: Lint = ', actirates%Lma_int_recrad(K), ' Lrec = ', actirates%Lma_recrad(K)
   Zintrecom = Zintrecom + actirates%Lma_int_recrad(K) 
   Zrecom = Zrecom + actirates%Lma_recrad(K)
   ! print*, 'Zrecom = ', Zrecom
   act_index = act_index + 1
  ELSE ! npoints = 0
   actirates%Lma_int_recrad(I) = 0.D0
   actirates%Lma_recrad(I) = 0.D0
  END IF ! npoints
 END DO
 ! STOP 'i_radion: testing'
END IF ! indexi > 1
! write(*,*) 'Zrecom = ', Zrecom
 write(*,*)  'photion_rates: Zion = ', Zion, ' Zrecom = ', Zrecom
! STOP 'i_radion: testing'
END SUBROUTINE i_radion
