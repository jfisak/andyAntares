! this subroutine find
  SUBROUTINE next_dyn_cell3(cell_numb,pack_index,next_cell)

  USE types
  IMPLICIT NONE

  ! input variables
  INTEGER                               :: cell_numb, pack_index
  ! output variables
  DOUBLE PRECISION, DIMENSION(3)        :: pos
  INTEGER                               :: next_cell
  ! parameter
  DOUBLE PRECISION                      :: t1, t2, t3, t4, t5, t6
  DOUBLE PRECISION                      :: tTest, tbound
  DOUBLE PRECISION                      :: solution, dist
  ! 
  DOUBLE PRECISION, DIMENSION(3)        :: corner, width, phot_pos, dir, testPos
  DOUBLE PRECISION, DIMENSION(3)        :: boundPos

  next_cell = 0
  ! we have to choose tTest small enough to be in next cell
  tTest = dist+minwidth/1.D66
  testPos = phot_pos + dir * tTest
  boundPos = phot_pos + dir * dist
  IF((testPos(1) <= -xmax) .OR. (testPos(1) >= xmax)) next_cell = -99
  IF((testPos(2) <= -ymax) .OR. (testPos(2) >= ymax)) next_cell = -99
  IF((testPos(3) <= -zmax) .OR. (testPos(3) >= zmax)) next_cell = -99
  IF(next_cell == 0) CALL find_dyn_cell(testPos,next_cell)
  ! now we will test which part of the cell were hit
  ! walls
!  print*, 'index = ', indx
  END SUBROUTINE next_dyn_cell3
