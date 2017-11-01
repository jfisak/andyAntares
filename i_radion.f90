SUBROUTINE i_radion(approx, indexe, indexi, leveli, current_mgi, act_pop, Zion, Zintrecom, Zrecom)!, actirates)
USE types
USE rates
IMPLICIT NONE

! input variables
INTEGER                                 :: approx, indexe, indexi, leveli
INTEGER                                 :: current_mgi, nrecom
DOUBLE PRECISION                        :: act_pop
! computing fields
INTEGER                                 :: npoints
INTEGER                                 :: I, Istart, J, Jstart, K
DOUBLE PRECISION, ALLOCATABLE           :: freq(:), cross(:), func(:)
DOUBLE PRECISION                        :: flux
DOUBLE PRECISION                        :: summ
DOUBLE PRECISION                        :: phot_cross
DOUBLE PRECISION                        :: freqt
DOUBLE PRECISION                        :: exci_energy, gr_exci_energy
DOUBLE PRECISION                        :: pop_number
DOUBLE PRECISION                        :: x
DOUBLE PRECISION                        :: flux_function
DOUBLE PRECISION                        :: stat_weight
! output variables
DOUBLE PRECISION                        :: Zion, Zrecom, Zintrecom
!CLASS(irates)                           :: actirates
INTEGER                                 :: OMP_GET_THREAD_NUM, my_rank

my_rank = OMP_GET_THREAD_NUM()


!write(*,*) 'i_radion: my_rank = ', my_rank, ' recrad = ', size(actirates%Lma_recrad), &
!           ' intrecrad = ', size(actirates%Lma_int_recrad), ' reccol = ', size(actirates%Lma_reccol), &
!           ' intreccol = ', size(actirates%Lma_int_reccol)
IF(indexi < elements(indexe)%atom_number) THEN
 npoints = SIZE(elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,:))
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
  freqt = (elements(indexe)%ions(indexi + 1)%levels(1)%exci_energy - &
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
  func(I) = cross(I) * flux / ( h * freq(I))
 ! print*, 'flux = ', flux, ' cross(I) = ', cross(I), ' func(I) = ', func(I), &
 !  ' h * freq = ', h * freq(I)
 END DO
 ! calculation of integral with the trapezoid rule
 DO I = Istart, npoints - 1
  summ = summ + (func(I) + func(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
  !print*, 'func(I) + func(I + 1) = ', func(I) + func(I + 1), &
  ! ' freq(I + 1) - freq(I) = ', freq(I + 1) - freq(I)
 END DO
 phot_cross = 4.D0 * pi * summ * act_pop
 !print*, 'phot_cross = ', phot_cross, ' summ = ', summ, ' act_pop = ', act_pop
 Zion = phot_cross * elements(indexe)%ions(indexi)%levels(leveli)%exci_energy
 DEALLOCATE(freq, cross, func)
ELSE
 Zion = 0.D0
END IF
 !print*, 'photion_rates: phot_cross = ', phot_cross
 !____________________________________________________________________________
 ! recombination
IF(indexi > 1) THEN
 Zintrecom = 0.D0
 Zrecom = 0.D0
 nrecom = SIZE(actirates%Lma_recrad)
 DO K = 1, nrecom
  npoints = SIZE(elements(indexe)%ions(indexi - 1)%levels(K)%photcros(1,:))
  IF (npoints /= 0) THEN
   ALLOCATE(freq(npoints), cross(npoints), func(npoints))
   !print*, 'photion_rates: indexe = ', indexe, ' indexi = ', indexi
   freq(1:npoints) = elements(indexe)%ions(indexi - 1)%levels(K)%photcros(1,1:npoints)
   cross(1:npoints) = elements(indexe)%ions(indexi - 1)%levels(K)%photcros(2,1:npoints)
   DO I = 1, npoints
    flux = flux_function(0,freq(I), T_eff)
    x = (h * freq(I)) / (BOLK * T_eff)
    !write(*,*) 'i_radion: f^3/c^2 = ', freq(I)**3/light_speed**2
    !write(*,*) 'i_radion: expr1 = ', cross(I)/(h * freq(I))! * ((2.0 * h * freq(I)**3.0)/light_speed**2.0 + flux) * exp(-x)
    !write(*,*) 'i_radion: expr2 = ', ((2.0 * h * freq(I)**3.0)/light_speed**2.0 + flux) * exp(-x)
    !write(*,*) 'i_radion: expr3 = ', ((2.0 * h * freq(I)**3.0)/light_speed**2.0 ) * exp(-x)
    IF(((2.0 * h * freq(I)**3.0) / light_speed**2.0 + flux) * exp(-x) > 1.D-150) THEN
     func(I) = cross(I) / (h * freq(I)) * &
      ((2.0 * h * freq(I)**3.0) / light_speed**2.0 + flux) * exp(-x)
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
    IF(freq(J) >= freqt) THEN
     Jstart = J
     EXIT
    IF(J == npoints) Jstart = 1
    END IF
   END DO
   !write(*,*) 'i_radion: Jstart = ', Jstart
   summ = 0
   DO J = Jstart, npoints - 1
    summ = summ + (func(J) + func(J + 1)) / 2.D0 * (freq(J + 1) - freq(J))
   END DO
   phot_cross = 4 * pi * act_pop * summ
   CALL populations(indexe, indexi - 1, K, current_mgi, pop_number)
   exci_energy = elements(indexe)%ions(indexi - 1)%levels(K)%exci_energy
   gr_exci_energy = elements(indexe)%ions(indexi)%levels(1)%exci_energy
   stat_weight = elements(indexe)%ions(indexi - 1)%levels(K)%stat_waight
   actirates%Lma_int_recrad(K) = pop_number * phot_cross * exci_energy * stat_weight
   actirates%Lma_recrad(K) = pop_number * phot_cross * stat_weight * (gr_exci_energy - exci_energy)
   !print*, 'photion_rates: actirates%Lma_recrad(K) = ', actirates%Lma_recrad(K), ' actirates%Lma_int_recrad = ', actirates%Lma_int_recrad(K)
   !actirates%Lma_recrad(K) = 0.D0
   Zintrecom = Zintrecom + actirates%Lma_int_recrad(K) 
   Zrecom = Zrecom + actirates%Lma_recrad(K)
   !print*, 'Zrecom = ', Zrecom
   DEALLOCATE(freq, cross, func)
  ELSE
   actirates%Lma_recrad(K) = 0.D0
   actirates%Lma_int_recrad(K) = 0.D0
  END IF
 END DO
ELSE
 Zrecom = 0.D0
 Zintrecom = 0.D0
END IF
 !print*, 'Zrecom = ', Zrecom
!print*, 'photion_rates: Zion = ', Zion, ' Zrecom = ', Zrecom

END SUBROUTINE i_radion
