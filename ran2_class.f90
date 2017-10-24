MODULE ran2_class

 IMPLICIT NONE

 INTEGER, PARAMETER               :: NTAB = 32

 TYPE, PUBLIC :: ranvar
 
  INTEGER, DIMENSION(NTAB)        :: iv
  INTEGER                         :: iy
  INTEGER                         :: idum2
  
  CONTAINS

  PROCEDURE, PUBLIC :: init => initialize_values
 END TYPE ranvar

 PRIVATE                          :: initialize_values

 CONTAINS

 SUBROUTINE initialize_values(this)
  CLASS(ranvar), POINTER:: this
  ! loop index
  INTEGER               :: ind
  
  DO ind = 1, NTAB
   this%iv(ind) = 0
  END DO
  this%iy = 0
  this%idum2 = 123456789
 END SUBROUTINE initialize_values


END MODULE ran2_class
