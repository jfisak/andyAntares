FUNCTION flux_function(approx, freq, temperature)
USE TYPES
IMPLICIT NONE

INTEGER                         :: approx
DOUBLE PRECISION                :: freq, temperature
DOUBLE PRECISION                :: flux_function


SELECT CASE(approx)
! the most stupid approximation: J is the Planck function
CASE(0)
 !write(*,*) 'expr = ', h * freq / (BOLK * temperature)
 IF( (h * freq / (BOLK * temperature)) < 7.D2) THEN
  flux_function = ( 2.D0 * h * freq**3 / light_speed**2 )  * &
   1.D0 / ( EXP( (h * freq / (BOLK * temperature) ) - 1.D0 ) )
 ELSE
  flux_function = 0.D0
 END IF
CASE DEFAULT
 STOP 'photion_rates: this approximation is not known'
END SELECT


END FUNCTION flux_function
