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
INTEGER                                 :: I, Istart, J, Jstart, K
DOUBLE PRECISION, ALLOCATABLE           :: freq(:), cross(:), func(:)
DOUBLE PRECISION                        :: gammaijk, photRate
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
! output variables
DOUBLE PRECISION                        :: Zion, Zrecom, Zintrecom
TYPE(irates)                           :: actirates
INTEGER                                 :: OMP_GET_THREAD_NUM, my_rank

my_rank = OMP_GET_THREAD_NUM()


!write(*,*) 'i_radion: my_rank = ', my_rank, ' recrad = ', size(actirates%Lma_recrad), &
!           ' intrecrad = ', size(actirates%Lma_int_recrad), ' reccol = ', size(actirates%Lma_reccol), &
!           ' intreccol = ', size(actirates%Lma_int_reccol)
el_dens = model_grid(current_mgi)%e_dens
temp = model_grid(current_mgi)%t
IF(indexi <= elements(indexe)%atom_number) THEN
 IF(ALLOCATED(elements(indexe)%ions(indexi)%levels(leveli)%photcros)) THEN
  npoints = SIZE(elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,:))
 ELSE
  npoints = 0
 END IF
 ! write(*,*) 'i_radion: npoints = ', npoints
ELSE
 npoints = 0
END IF
!print*, 'photion_rates: npoints = ', npoints
! there are no data for photoionization cross section available
! the rates are equal to zero
!print*, 'photion_rates: npoints = ', npoints
IF(npoints /= 0) THEN
 ALLOCATE(freq(npoints), cross(npoints), func(npoints))
 T_eff = model_grid(current_mgi)%T
 freq(1:npoints) = elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,1:npoints)
 cross(1:npoints) = elements(indexe)%ions(indexi)%levels(leveli)%photcros(2,1:npoints)
 freqt = (MINVAL(elements(indexe)%ions(indexi + 1)%levels(:)%exci_energy) - &
  elements(indexe)%ions(indexi)%levels(leveli)%exci_energy) / h
 ! looking for starting point
 DO I = 1, npoints
  IF(freq(I) >= freqt) THEN
   Istart = I
   EXIT
  END IF
 END DO
 summ = 0.D0
 DO I = Istart, npoints
   flux = flux_function(0,freq(I), T_eff)
  ! func(I) = cross(I) * flux / ( h * freq(I))
  func(I) = flux * cross(I) / (h * freq(I)) * (1 - exp(-(h * freq(I) / (BOLK * temp))))
 ! print*, 'flux = ', flux, ' cross(I) = ', cross(I), ' func(I) = ', func(I), &
 !  ' h * freq = ', h * freq(I)
 END DO
 ! calculation of integral with the trapezoid rule
 DO I = Istart, npoints - 1
  ! summ = summ + (func(I) + func(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
  summ = summ + (func(I) + func(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
  !print*, 'func(I) + func(I + 1) = ', func(I) + func(I + 1), &
  ! ' freq(I + 1) - freq(I) = ', freq(I + 1) - freq(I)
 END DO
 CALL saha_factor(indexe, indexi + 1, leveli, current_mgi, el_dens, sfactor)
 CALL populations(indexe, indexi + 1, 1, current_mgi, up_pop)
 gammaijk = 4.D0 * pi * summ
 phot_cross = 4.D0 * pi * summ * act_pop * sfactor
 photRate = act_pop * gammaijk
  ! write(*,*) 'i_radion: phot_cross = ', phot_cross, ' summ = ', summ, ' act_pop = ', act_pop, ' photRate = ', photRate
  ! write(*,*) 'i_radion: up_pop = ', up_pop
  ! write(*,*) 'exci_energy = ', elements(indexe)%ions(indexi)%levels(leveli)%exci_energy
 IF(photRate < 0.D0) STOP 'i_radion: photRate < 0'
 Zion = photRate * elements(indexe)%ions(indexi)%levels(leveli)%exci_energy
 ! write(*,*) 'i_radion: Zion = ', Zion
 DEALLOCATE(freq, cross, func)
ELSE
 Zion = 0.D0
END IF
 ! write(*,*) 'i_radion: Zion = ', Zion
 ! print*, 'photion_rates: phot_cross = ', phot_cross
 !____________________________________________________________________________
 ! recombination
IF(indexi > 1) THEN
 Zintrecom = 0.D0
 Zrecom = 0.D0
 nrecom = SIZE(actirates%Lma_recrad)
 ! write(*,*) 'i_radion: nrecom = ', nrecom
 DO K = 1, nrecom
  IF(nrecom == 1) EXIT
   IF(ALLOCATED(elements(indexe)%ions(indexi - 1)%levels(K)%photcros)) THEN
    npoints = SIZE(elements(indexe)%ions(indexi - 1)%levels(K)%photcros(1,:))
   ELSE
    npoints = 0
   END IF
  IF (npoints /= 0) THEN
   ALLOCATE(freq(npoints), cross(npoints), func(npoints))
   !print*, 'photion_rates: indexe = ', indexe, ' indexi = ', indexi
   freq(1:npoints) = elements(indexe)%ions(indexi - 1)%levels(K)%photcros(1,1:npoints)
   cross(1:npoints) = elements(indexe)%ions(indexi - 1)%levels(K)%photcros(2,1:npoints)
   DO I = 1, npoints
    ! flux = flux_function(0,freq(I), T_eff)
    x = (h * freq(I)) / (BOLK * T_eff)
    IF(((2.0 * h * freq(I)**3.0) / light_speed**2.0 + flux) * exp(-x) > 1.D-150) THEN
     func(I) = cross(I) / (h * freq(I)) * &
      (2.0 * h * freq(I)**3.0) / light_speed**2.0 * exp(-x)
    ELSE
     func(I) = 0.D0
    END IF
   END DO
   freqt = (elements(indexe)%ions(indexi)%levels(leveli)%exci_energy - &
    elements(indexe)%ions(indexi - 1)%levels(K)%exci_energy) / h
   !write(*,*) 'i_radion: excienergy1 = ', &
   ! elements(indexe)%ions(indexi)%levels(leveli)%exci_energy, &
   ! ' excienergy2 =  ',  &
   ! elements(indexe)%ions(indexi - 1)%levels(K)%exci_energy, &
   ! 'freq = ', freqt
   ! looking for starting point
   DO J = 1, npoints
    ! write(*,*) 'i_radion: J = ', J
    IF(freq(J) >= freqt) THEN
     Jstart = J
     EXIT
    END IF
    IF(J == npoints) Jstart = 1
   END DO
   summ = 0
   DO J = Jstart, npoints - 1
    summ = summ + (func(J) + func(J + 1)) / 2.D0 * (freq(J + 1) - freq(J))
   END DO
   phot_cross = 4.0 * pi * summ
   ! write(*,*) 'i_radion: calling populations...'
   CALL populations(indexe, indexi, 1, current_mgi, pop_number)
   CALL saha_factor(indexe, indexi, K, current_mgi, el_dens, sfactor)
   ! CALL saha_boltzmann_factor(indexe, indexi, temp, sb_factor, too_large
   exci_energy = elements(indexe)%ions(indexi - 1)%levels(K)%exci_energy
   gr_exci_energy = MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
   actirates%Lma_int_recrad(K) = pop_number * phot_cross * sfactor * exci_energy
   actirates%Lma_recrad(K) = pop_number * phot_cross * sfactor *&
    (gr_exci_energy - exci_energy)
   ! write(*,*) 'i_radion: indexe = ', indexe, ' indexi - 1 = ', indexi - 1, 'indexl = ', K
   ! write(*,*) 'i_radion: exci_energy = ', exci_energy, ' gr_exci_energy = ', gr_exci_energy
   ! write(*,*) 'i_radion: pop_number = ', pop_number, ' phot_cross = ', phot_cross
   ! write(*,*) 'i_radion: Lma_int_recrad(', K, ') = ', actirates%Lma_int_recrad(K)
   Zintrecom = Zintrecom + actirates%Lma_int_recrad(K) 
   Zrecom = Zrecom + actirates%Lma_recrad(K)
   ! print*, 'Zrecom = ', Zrecom
   DEALLOCATE(freq, cross, func)
  ELSE
   actirates%Lma_recrad(K) = 0.D0
   actirates%Lma_int_recrad(K) = 0.D0
  END IF
 END DO
ELSE ! indexi > 1
 Zrecom = 0.D0
 Zintrecom = 0.D0
END IF ! indexi > 1
! print*, 'Zrecom = ', Zrecom
!print*, 'photion_rates: Zion = ', Zion, ' Zrecom = ', Zrecom

END SUBROUTINE i_radion
