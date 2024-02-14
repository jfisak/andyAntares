MODULE ran2_class

 IMPLICIT NONE

 INTEGER, PARAMETER               :: NTAB = 32
 INTEGER, PARAMETER               :: IM1 = 2147483563, &
                                     IM2 = 2147483399, &
                                     IMM1 = IM1 - 1,   &
                                     IA1 = 40014, &
                                     IA2 = 40692, &
                                     IQ1 = 53668, &
                                     IQ2 = 52774, &
                                     IR1 = 12211, &
                                     IR2 = 3791, &
                                     NDIV = 1 + INT(IMM1 / NTAB)
 DOUBLE PRECISION, PARAMETER      :: AM = 1. / DBLE(IM1), &
                                     EPS = 1.2D-7, &
                                     RNMX = 1. - EPS

 TYPE, PUBLIC :: ranvar
  INTEGER, DIMENSION(NTAB)        :: iv
  INTEGER                         :: iy
  INTEGER                         :: idum2
  
  !CONTAINS

  !PROCEDURE, PUBLIC :: init_val => initialize_values
 END TYPE ranvar

 ! global member of the ranvar class
 CLASS(ranvar), POINTER      :: ranum 
 !PRIVATE                          :: initialize_values
 INTERFACE ranvar
  MODULE PROCEDURE initialize_values
 END INTERFACE

 CONTAINS

 FUNCTION initialize_values()
  IMPLICIT NONE
  CLASS(ranvar), POINTER:: initialize_values
  ! loop index
  INTEGER               :: ind
  
  DO ind = 1, NTAB
   initialize_values%iv(ind) = 0
  END DO
  initialize_values%iy = 0
  initialize_values%idum2 = 123456789
 END FUNCTION initialize_values



END MODULE ran2_class
