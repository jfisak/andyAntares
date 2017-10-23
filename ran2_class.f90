MODULE ran2_class

 IMPLICIT NONE


 TYPE, PRIVATE :: ranvar
 
  INTEGER                         :: NTAB
  INTEGER, ALLOCATABLE            :: iv(:)
  INTEGER                         :: iy
  INTEGER                         :: idum2
  
  CONTAINS

  PROCEDURE, PUBLIC :: init => initialize_values
 END TYPE ranvar

 PRIVATE                          :: initialize_values

 CONTAINS

 SUBROUTINE initialize_values(this)
  CLASS(ranvar)         :: this
  ! loop index
  INTEGER               :: ind
  INTEGER               :: NTAB

  NTAB = 32
  this%NTAB = NTAB
  IF(.NOT. ALLOCATED(this%iv)) ALLOCATE(this%iv(NTAB))
  DO ind = 1, NTAB
   this%iv(ind) = 0
  END DO
  this%iy = 0
  this%idum2 = 123456789
 END SUBROUTINE initialize_values

END MODULE ran2_class
