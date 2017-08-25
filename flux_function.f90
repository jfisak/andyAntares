FUNCTION flux_function(approx, freq, temperature)
USE TYPES
IMPLICIT NONE

INTEGER                         :: approx
DOUBLE PRECISION                :: freq, temperature
DOUBLE PRECISION                :: flux_function


SELECT CASE(approx)
! the most stupid approximation: J is the Planck function
CASE(0)
 flux_function = ( 2.D0 * h * freq**3 / light_speed**2 )  * &
  1.D0 / ( EXP( (h * freq / (BOLK * temperature) ) - 1.D0 ) )
CASE DEFAULT
 STOP 'photion_rates: this approximation is not known'
END SELECT


END FUNCTION flux_function
