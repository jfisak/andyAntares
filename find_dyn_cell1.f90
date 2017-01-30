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
! parameters of subcells of dyngrid ijk
DOUBLE PRECISION, DIMENSION(3)                  :: subcells_width
INTEGER                                         :: subind_x, subind_y, subind_z


! firstly we can compute which basic cell this point contains
bcell(1) = FLOOR(pos(1)/basic_cell_width(1) + dble(nx_cell)/2.D0) + 1
bcell(2) = FLOOR(pos(2)/basic_cell_width(2) + dble(ny_cell)/2.D0) + 1
bcell(3) = FLOOR(pos(3)/basic_cell_width(3) + dble(nz_cell)/2.D0) + 1
! index of the given basic cell
bindex = (bcell(1) - 1) * ny_cell * nz_cell + (bcell(2) - 1) * nz_cell + bcell(3)
! initial setting of the local variable corresponding to the actual cell
actCell = bindex
!print*, 'find_dyn_cell1: actCell = ', actCell
! if there is no dynamical cell in the given basic cell

IF(dyn_cell(actCell)%up_cell == 0) THEN
 actual_cell = actCell
 RETURN 
! we will move upper in the dyncell tree otherwise
END IF
SELECT CASE(dyngrid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dyngrid 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
! we are looking for the given cell in dyncell tree
actCell = dyn_cell(actCell)%up_cell
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dyngrid ijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 subcells_width = dyn_cell(dyn_cell(actCell)%up_cell)%width
 subind_x = FLOOR((dyn_cell(actcell)%corner(1) - pos(1))/subcells_width(1))
 subind_y = FLOOR((dyn_cell(actcell)%corner(2) - pos(2))/subcells_width(2))
 subind_z = FLOOR((dyn_cell(actcell)%corner(3) - pos(3))/subcells_width(3))
 print*, 'next_cell_up:', subind_x, subind_y, subind_z
 actual_cell = dyn_cell(actcell)%up_cell + &
        subcells_width(1) * subcells_width(2) * (subind_x - 1) + &
        subcells_width(1) * (subind_y - 1) + subind_z - 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! default case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE DEFAULT
 STOP 'find_dyn_cell1: choice of the dyngrid type is not known'
END SELECT

END SUBROUTINE find_dyn_cell1
