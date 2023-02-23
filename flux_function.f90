FUNCTION flux_function(approx, freq, temperature, cur_r)
USE TYPES
USE constants
IMPLICIT NONE

INTEGER                         :: approx
DOUBLE PRECISION                :: freq, temperature
DOUBLE PRECISION                :: flux_function
DOUBLE PRECISION                :: W, cur_r


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
! diluted blackbody
CASE(1)
 ! the dilution factor
 W = 5.D-1 * (1.0-sqrt(1.0-(R_star/cur_r)**2.0))
 ! temperature = 1.E4
 IF( (h * freq / (BOLK * temperature)) < 7.D2) THEN
  flux_function = W * ( 2.D0 * h * freq**3 / light_speed**2 )  * &
   1.D0 / ( EXP( (h * freq / (BOLK * temperature) ) - 1.D0 ) )
 ELSE
  flux_function = 0.D0
 END IF
 
CASE DEFAULT
 STOP 'photion_rates: this approximation is not known'
END SELECT


END FUNCTION flux_function
