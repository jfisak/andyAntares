! this subroutine find
  SUBROUTINE next_dyn_cell(pack_index,t,tTest,next_cell)

  USE types
  IMPLICIT NONE

  ! input variables
  INTEGER                               :: cell_numb, pack_index
  ! output variables
  DOUBLE PRECISION, DIMENSION(3)        :: pos
  INTEGER                               :: next_cell
  ! parameter
  DOUBLE PRECISION                      :: t
  DOUBLE PRECISION                      :: tTest, tbound
  DOUBLE PRECISION                      :: solution, dist
  ! 
  DOUBLE PRECISION, DIMENSION(3)        :: corner, width, phot_pos, dir, testPos
  DOUBLE PRECISION, DIMENSION(3)        :: boundPos

  boundPos = package(pack_index)%pos + package(pack_index)%dir * t
  IF(boundPos(1) == -xmax * r_sun .OR. boundPos(1) == xmax * r_sun) next_cell = -99
  IF(boundPos(2) == -ymax * r_sun .OR. boundPos(2) == ymax * r_sun) next_cell = -99
  IF(boundPos(3) == -zmax * r_sun .OR. boundPos(3) == zmax * r_sun) next_cell = -99
  testPos = package(pack_index)%pos + package(pack_index)%dir * tTest
  print*, 'next_cell = ', next_cell
  IF(next_cell /= -99) CALL find_dyn_cell1(testPos,next_cell)
  END SUBROUTINE next_dyn_cell
