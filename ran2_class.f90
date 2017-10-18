MODULE ran2_class

 IMPLICIT NONE


 TYPE, PRIVATE :: ranvar
 
  INTEGER, PARAMETER              :: NTAB = 32
  INTEGER, DIMENSION(NTAB)        :: iv
  INTEGER                         :: iy
  INTEGER                         :: idum2
  
  CONTAINS

  PROCEDURE, PUBLIC :: initialize_values
 END TYPE ranvar


 SUBROUTINE initialize_values()
  CLASS(ranvar)         :: this
  ! loop index
  INTEGER               :: ind

  DO ind = 1, NTAB
   this%iv(ind) = 0
  END DO
  this%iy = 0
  this%idum2 = 123456789
 END SUBROUTINE

END MODULE ran2_class
