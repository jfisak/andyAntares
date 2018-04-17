SUBROUTINE col_deexcitations_event(cur_element, cur_ion, line)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: cur_element, cur_ion, line
INTEGER                         :: nlns
INTEGER, ALLOCATABLE            :: linecoltransitions(:), Lcoldeex(:)

! firstly we have to compute number of possible collisional transitions from last_line
nlns = 0
rand = ran2(idum)
DO I = 1, ntransitions
 ! we are interested only in the transitions for the given atom
 IF(linelist(I)%indexe == element_index .AND. linelist(I)%indexi == ion_index) THEN
  ! transitions to a lower level
  IF(linelist(I)%upper == line) THEN
   nlns = nlns + 1
  END IF
 END IF
END DO
ALLOCATE(linecoltransitions(nlns), Lcoldeex(nlns))
Zcoldeexc = 0.D0
J = 0
DO I = 1, ntransitions
 IF(linelist(I)%indexe == element_index .AND. linelist(I)%indexi == ion_index) THEN
  ! transitions to a lower level
  IF(linelist(I)%upper == line) THEN
   J = J + 1
   linecoltransitions(J) = I
   actVal = linelist(I)%A_ul
   Lcoldeex(J) = actVal
   Zcoldeexc = Zcoldeexc + actVal
  END IF
 END IF
END DO
! we have to now choose the given process of deexcitation
rand = rand * Zcoldeexc




END SUBROUTINE
