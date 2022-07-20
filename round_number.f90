! a functinon for a simple calculation of a round number
FUNCTION round_number(input_number)

IMPLICIT NONE

DOUBLE PRECISION                       :: input_number, round_number, zbytek, vysledek

zbytek = MODULO(input_number, 1.0)

IF( zbytek == 0.0) THEN
 vysledek = input_number
ELSE IF (zbytek < 0.5) THEN
 vysledek = FLOOR(input_number)
ELSE IF (zbytek >= 0.5) THEN
 vysledek = CEILING(input_number)
END IF

 round_number = vysledek
 RETURN

END FUNCTION
