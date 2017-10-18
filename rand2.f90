MODULE ran2_class

 IMPLICIT NONE


 TYPE, PUBLIC :: ranvar
 
  INTEGER, PARAMETER              :: NTAB = 32
  INTEGER, DIMENSION(32)          :: iv = NTAB * 0
  INTEGER                         :: iy = 0
  INTEGER                         :: idum2=123456789
 END TYPE ranvar


END MODULE ran2_class
