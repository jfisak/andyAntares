! SBR will find an element index for the given atomic number
SUBROUTINE find_element_index(Z, indexe)
USE types
IMPLICIT NONE

INTEGER                                 :: Z
INTEGER                                 :: indexe
INTEGER                                 :: I

DO I = 1, n_elements
 IF(Z == elements(I)%atom_number) THEN
  indexe = I
  RETURN
 END IF
END DO

STOP 'find_element_index: no valid element was found'

END SUBROUTINE find_element_index
