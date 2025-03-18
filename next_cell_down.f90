! finds the next cell in the lower levels of the propGrid tree
!
! INPUT: pack_index(INT): packet index
! OUTPUT: next_cell(INT): index of the next cell
!
SUBROUTINE next_cell_down(pack_index, next_cell)

USE types
USE constants
USE dummypacket

IMPLICIT NONE

! input variables
INTEGER                                 :: pack_index, dummypack_index
! outpu variables
INTEGER                                 :: next_cell
! local variable next cell
INTEGER                                 :: n_cell
! cross the surface
INTEGER                                 :: cross
! actual cell
INTEGER                                 :: act_cell

DOUBLE PRECISION                        :: calc_dist
INTEGER                                 :: cell_numb
INTEGER                                 :: n_pos, n_neg, n_zer, n_par


IF(pack_index > SIZE(package)) THEN
 dummypack_index = pack_index - SIZE(package)
 act_cell = dummypackage(dummypack_index)%cell_numb
 cross = dummypackage(dummypack_index)%next_cross
ELSE
 act_cell = package(pack_index)%cell_numb
 cross = package(pack_index)%next_cross
END IF

IF(cross <= 0 ) THEN
 ! cell_numb and calc_dist are not important in this case
 CALL bound_dist(pack_index, cell_numb, calc_dist, n_pos, n_neg, n_zer, n_par)
 ! cross = package(pack_index)%next_cross
 IF(cross < 0) STOP 'next_cell_down: next cross is impossible to find'
END IF

DO
 n_cell = dyn_cell(act_cell)%neighbor(cross)
 IF(n_cell == 0) THEN
  IF(dyn_cell(act_cell)%down_cell == 0) STOP 'next_cell: no cell was found'
  act_cell = dyn_cell(act_cell)%down_cell
 ELSE 
  next_cell = n_cell
  IF(next_cell > SIZE(dyn_cell)) THEN
   write(*,*) 'next_cell_down: next_cell = ', next_cell, ' > number of propGrid cells = ', SIZE(dyn_cell)
   CALL abort()
  END IF
  EXIT
 END IF
END DO

END SUBROUTINE next_cell_down
