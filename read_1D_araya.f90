SUBROUTINE read_1D_araya()

! file format
! Teff
! R_star
! R/R_* -- V/km/s -- -- rho/g/cm^3 -- -- -- -- -- -- -- -- -- -- --

USE types
IMPLICIT NONE

INTEGER                                 :: n_mg_points

INTEGER                                 :: cur_line, reading_grid

DOUBLE PRECISION                        :: effective_temperature, stellar_radius
DOUBLE PRECISION                        :: radius, density, djunk
DOUBLE PRECISION                        :: velocity

CHARACTER(LEN=100)                      :: junk

add_mg = 1

n_mg_points = 0
OPEN(38, FILE=inputmodelFile)

 READ(38,*) junk
 READ(38,*) junk
 ! number of model points
 DO 
  READ(38,*,iostat=reading_grid) junk
  IF(reading_grid /= 0) EXIT
  n_mg_points = n_mg_points + 1
 END DO

 ALLOCATE(model_grid(n_mg_points + add_mg))

 REWIND(38)

 READ(38,*) effective_temperature
 READ(38,*) stellar_radius


 DO cur_line = 1, n_mg_points
  READ(38,*) radius, djunk, velocity, djunk, djunk, density, junk
  model_grid
  write(*,*) 'read_1D_araya: r = ', radius, ' v = ', velocity, ' density = ', density
 END DO
 STOP 'read_1D_araya: testing'









CLOSE(38)

END SUBROUTINE read_1D_araya
