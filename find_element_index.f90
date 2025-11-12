! SBR will find an element index for the given atomic number
!
! INPUT: atomic_index(INT) -- atomic index
!        indexe(INT) -- element index
! OUTPUT: NONE
!
! 1x RETURN point
!
SUBROUTINE find_element_index(atomic_index, indexe)
USE types
USE constants
IMPLICIT NONE

INTEGER                                 :: atomic_index
INTEGER                                 :: indexe
INTEGER                                 :: ind_I

DO ind_I = 1, n_elements
 IF(atomic_index == elements(ind_I)%atom_number) THEN
  indexe = ind_I
  ! RETURN point
  RETURN
 END IF
END DO

STOP 'find_element_index: no valid element was found'

END SUBROUTINE find_element_index
