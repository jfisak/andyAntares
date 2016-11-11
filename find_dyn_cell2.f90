! this subroutine will find all dynamic cells
! containing the given point (this is the main
! difference between this sbr and sbr find_dyn_cell1
SUBROUTINE find_dyn_cell2(pos, act_dyn_cell, basic_cell, next_cell)

USE types
IMPLICIT NONE

! input variables
! position of the given point
DOUBLE PRECISION, DIMENSION(3)                  :: pos
! index of the next basic cell
INTEGER                                         :: basic_cell
! actual dynamic cell of the photon
INTEGER                                         :: act_dyn_cell
! output variables
INTEGER                                         :: next_cell
! basic cell variables
        INTEGER, DIMENSION(3)                   :: bcell
        INTEGER                                 :: bindex
! local actual cell
        INTEGER                                 :: actCell
! cells we are looking for
        INTEGER                                 :: foundCells
        INTEGER, DIMENSION(8)                   :: cells

actCell = basic_cell
! there is no dynamical cell we are looking for
! so the next cell we are looking for is the
! basic cell with index basic_cell
IF(dyn_cell(actCell)%up_cell == 0) THEN
 next_cell = actCell    
 RETURN                 
ELSE
 actCell = dyn_cell(actCell)%up_cell
END IF
!_____________________________________________________________________
! we are looking for the given cells in dyncell tree !!!!!!!!!!!!!!!!!
!_____________________________________________________________________
DO
 ! did we found the given cell containing the given point?
 IF((pos(1) .GE. dyn_cell(actCell)%corner(1)) .AND. &
    (pos(1) .LE. dyn_cell(actCell)%corner(1) + dyn_cell(actCell)%width(1)) .AND. &
    (pos(2) .GE. dyn_cell(actCell)%corner(2)) .AND. &
    (pos(2) .LE. dyn_cell(actCell)%corner(2) + dyn_cell(actCell)%width(2)) .AND. &
    (pos(3) .GE. dyn_cell(actCell)%corner(3)) .AND. &
    (pos(3) .LE. dyn_cell(actCell)%corner(3) + dyn_cell(actCell)%width(3))) THEN
  ! we have found a cell containing the given point
  ! is this cell on the top of the dyncell tree?
  IF(dyn_cell(actCell)%up_cell == 0) THEN
   foundCells = foundCells + 1
   cells(foundCells) = actCell
  ! we have to move to the higher level of the dyncell tree
  ELSE
   actCell = dyn_cell(actCell)%up_cell
  END IF
 ! we are not in the right cell, we have to move to the next cell
 ! in the dynamical cells tree
 END IF
 IF(actCell + 1 .LE. dyn_cell(dyn_cell(actCell)%down_cell)%up_cell + 7) THEN
  actCell = actCell + 1
 ELSE
  ! we move to the lower level of the dyncell grid
  actCell = dyn_cell(actCell)%down_cell
 END IF
 IF(actCell == basic_cell) EXIT
END DO
!_____________________________________________________________________
!_____________________________________________________________________
!_____________________________________________________________________
! now we have a set of possible cells we have to choose which one is the right cell
print*, 'find_dyn_cell2: foundCells = ', foundCells
IF(foundCells == 0) THEN
 STOP 'no cell was found'
ELSE IF(foundCells == 1) THEN
 next_cell = cells(1)
ELSE IF(foundCells == 2) THEN
 if (cells(1) == act_dyn_cell) then
  next_cell = cells(2)
 else if(cells(2) == act_dyn_cell) then
  next_cell = cells(1)
 end if
print*, 'find_dyn_cell2: next_cell = ', next_cell
! we don't have a solution for more cells now
! will be added in the future
ELSE
 STOP 'this number of cells is not now known'
END IF

END SUBROUTINE find_dyn_cell2
