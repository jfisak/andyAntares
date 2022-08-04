! this sbr analyses input from the file input.dat
SUBROUTINE analyse_input()

USE types
IMPLICIT NONE



IF (dyngrid /= 0) CALL virtual_particles(model_type)
 xmax = R_inf + 0.5 * R_sun
 ymax = R_inf + 0.5 * R_sun 
IF(model_type == 1) THEN
 zmax = R_inf + 0.5 * R_sun
ELSE IF (model_type == 2 .AND. inputmodel == 1) THEN
 zmax = Z_inf! + R_sun
ELSE IF (model_type == 2 .AND. inputmodel == 2) THEN
 zmax = R_inf + 0.5 * R_sun
ELSE
 STOP 'main: non-known model type'
END IF




END SUBROUTINE analyse_input
