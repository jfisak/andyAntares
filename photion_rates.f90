SUBROUTINE photion_rates(approx, indexe, indexi, leveli, current_mgi, act_pop, Zcross)
USE types
IMPLICIT NONE

! input variables
INTEGER                                 :: approx, indexe, indexi, leveli
INTEGER                                 :: current_mgi
DOUBLE PRECISION                        :: act_pop
! computing fields
INTEGER                                 :: npoints
INTEGER                                 :: I
DOUBLE PRECISION, ALLOCATABLE           :: freq(:), cross(:), func(:)
DOUBLE PRECISION                        :: planck
DOUBLE PRECISION                        :: summ
DOUBLE PRECISION                        :: phot_cross
! output variables
DOUBLE PRECISION                        :: Zcross

SELECT CASE(approx)
! the most stupid approximation: J is the Planck function
CASE(0)
 npoints = SIZE(elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,:))
 !print*, 'photion_rates: npoints = ', npoints
 IF(npoints == 0) THEN
  phot_cross = 0.D0
  RETURN
 END IF
 ALLOCATE(freq(npoints), cross(npoints), func(npoints))
 T_eff = model_grid(current_mgi)%T
 freq(1:npoints) = elements(indexe)%ions(indexi)%levels(leveli)%photcros(1,1:npoints)
 cross(1:npoints) = elements(indexe)%ions(indexi)%levels(leveli)%photcros(2,1:npoints)
 summ = 0.D0
 DO I = 1, npoints
  planck = ( 2.D0 * h * freq(I)**3 / light_speed**2 )  * &
          1.D0 / ( EXP( (h * freq(I) / (BOLK * T_eff) ) - 1.D0 ) )
  func(I) = cross(I) * planck / ( h * freq(I))
 END DO
 DO I = 1, npoints - 1
  summ = summ + (func(I) + func(I + 1)) / 2.D0 * (freq(I + 1) - freq(I))
 END DO
 phot_cross = 4 * pi * summ * act_pop
 Zcross = phot_cross * elements(indexe)%ions(indexi)%levels(leveli)%exci_energy
 print*, 'photion_rates: phot_cross = ', phot_cross
CASE DEFAULT
 STOP 'photion_rates: this approximation is not known'
END SELECT

END SUBROUTINE photion_rates
