! this subroutine finds a dynamic cell, where is the given
! photon located
! this subroutine expects a position out of bounds of the dyncells
SUBROUTINE find_dyn_cell1(pos,actual_cell)

USE types
IMPLICIT NONE

! input variables
        DOUBLE PRECISION, DIMENSION(3)          :: pos
! subroutine returns number of the actual cell
        INTEGER                                 :: actual_cell
! basic cell variables
        INTEGER, DIMENSION(3)                   :: bcell
        INTEGER                                 :: bindex
! local actual cell
        INTEGER                                 :: actCell

! firstly we can compute which basic cell this point contains
bcell(1) = FLOOR(pos(1)/basic_cell_width(1) + dble(nx_cell)/2.D0) + 1
bcell(2) = FLOOR(pos(2)/basic_cell_width(2) + dble(ny_cell)/2.D0) + 1
bcell(3) = FLOOR(pos(3)/basic_cell_width(3) + dble(nz_cell)/2.D0) + 1
! index of the given basic cell
bindex = (bcell(1) - 1) * ny_cell * nz_cell + (bcell(2) - 1) * nz_cell + bcell(3)
! initial setting of the local variable corresponding to the actual cell
actCell = bindex
print*, 'find_dyn_cell1: actCell = ', actCell
! if there is no dynamical cell in the given basic cell
IF(dyn_cell(actCell)%up_cell == 0) THEN
 actual_cell = actCell
 RETURN 
! we will move upper in the dyncell tree otherwise
ELSE
 actCell = dyn_cell(actCell)%up_cell
END IF
! we are looking for the given cell in dyncell tree
DO
 ! did we found the given cell containing the given point?
 IF((pos(1) .GE. dyn_cell(actCell)%corner(1)) .AND. (pos(1) .LE. dyn_cell(actCell)%corner(1) + dyn_cell(actCell)%width(1)) .AND. &
    (pos(2) .GE. dyn_cell(actCell)%corner(2)) .AND. (pos(2) .LE. dyn_cell(actCell)%corner(2) + dyn_cell(actCell)%width(2)) .AND. &
    (pos(3) .GE. dyn_cell(actCell)%corner(3)) .AND. (pos(3) .LE. dyn_cell(actCell)%corner(3) + dyn_cell(actCell)%width(3))) THEN
  ! we have found a cell containing the given point
  ! is this cell on the top of the dyncell tree?
  IF(dyn_cell(actCell)%up_cell == 0) THEN
   actual_cell = actCell
   EXIT
  ! we have to move to the higher level of the dyncell tree
  ELSE
   actCell = dyn_cell(actCell)%up_cell
  END IF
 ! we are not in the right cell, we have to move to the next cell
 ! in the dynamical cells tree
 ELSE
  IF(actCell + 1 .LE. dyn_cell(dyn_cell(actCell)%down_cell)%up_cell + 7) THEN
   actCell = actCell + 1
  ELSE
   STOP 'error in the next dynamical cell calculating'
  END IF
 END IF
END DO

END SUBROUTINE find_dyn_cell1
