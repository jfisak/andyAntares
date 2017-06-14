SUBROUTINE photion_rates(approx, indexe, indexi, leveli, current_mgi, act_pop, Zion, nrecom, Lrecom, Zrecom)
USE types
IMPLICIT NONE

! input variables
INTEGER                                 :: approx, indexe, indexi, leveli
INTEGER                                 :: current_mgi, nrecom
DOUBLE PRECISION, DIMENSION(nrecom)     :: Lrecom
DOUBLE PRECISION                        :: act_pop
! computing fields
INTEGER                                 :: npoints
INTEGER                                 :: I, Istart
DOUBLE PRECISION, ALLOCATABLE           :: freq(:), cross(:), func(:)
DOUBLE PRECISION                        :: flux
DOUBLE PRECISION                        :: summ
DOUBLE PRECISION                        :: phot_cross, coll_cross
DOUBLE PRECISION                        :: freqt
DOUBLE PRECISION                        :: exci_energy
DOUBLE PRECISION                        :: pop_number
DOUBLE PRECISION                        :: x
DOUBLE PRECISION                        :: flux_function
DOUBLE PRECISION                        :: stat_weight
! output variables
DOUBLE PRECISION                        :: Zion, Zrecom

npoints = SIZE(elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,:))
!print*, 'photion_rates: npoints = ', npoints
! there are no data for photoionization cross section available
! the rates are equal to zero
IF(npoints == 0) THEN
 print*, 'no valid data for phion cs calculation...'
 Zion = 0.D0
 Zrecom = 0.D0
 RETURN
END IF
ALLOCATE(freq(npoints), cross(npoints), func(npoints))
T_eff = model_grid(current_mgi)%T
freq(1:npoints) = elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,1:npoints)
cross(1:npoints) = elements(indexe)%ions(indexi)%levels(leveli)%photcros(2,1:npoints)
IF(indexi < elements(indexe)%atom_number) THEN
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
ELSE
 Zion = 0.D0
END IF
 !print*, 'photion_rates: phot_cross = ', phot_cross
 !____________________________________________________________________________
 ! recombination
IF(indexi > 1) THEN
 DO I = 1, npoints
  flux = flux_function(0,freq(I), T_eff)
  x = (h * freq(I)) / (BOLK * T_eff)
  func(I) = cross(I) / (h * freq(I)) * &
   ((2 * h * freq(I)**3.0) / light_speed**2.0 + flux) * exp(-x)
  !print*, 'func(I) = ', func(I)
 END DO
 summ = 0
 DO I = 1, npoints - 1
  summ = summ + (func(I) + func(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
 END DO
 coll_cross = 4 * pi * act_pop * summ
 Zrecom = 0.D0
 DO I = 1, nrecom
  CALL populations(indexe, indexi - 1, I, current_mgi, pop_number)
  exci_energy = elements(indexe)%ions(indexi - 1)%levels(I)%exci_energy
  stat_weight = elements(indexe)%ions(indexi - 1)%levels(I)%stat_waight
  Lrecom(I) = pop_number * coll_cross * exci_energy * stat_weight
  Lrecom(I) = 0.D0
  Zrecom = Zrecom + Lrecom(I) 
 END DO
ELSE
 Zrecom = 0.D0
END IF
 !print*, 'Zrecom = ', Zrecom
print*, 'photion_rates: Zion = ', Zion

END SUBROUTINE photion_rates
