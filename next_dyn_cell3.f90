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

  

















  END SUBROUTINE next_dyn_cell3
