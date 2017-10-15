MODULE rand2



 IMPLICIT NONE
 
 !$THREATPRIVATE(NTAB, iv,iy,idum2)
 INTEGER, PARAMETER              :: NTAB = 32
 INTEGER, DIMENSION(32)          :: iv = NTAB * 0
 INTEGER                         :: iy = 0
 INTEGER                         :: idum2=123456789


END MODULE rand2
