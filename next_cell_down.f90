SUBROUTINE next_cell_down(pack_index, next_cell)

USE types

IMPLICIT NONE

! input variables
INTEGER                                 :: pack_index
! outpu variables
INTEGER                                 :: next_cell
! local variable next cell
INTEGER                                 :: n_cell
! cross the surface
INTEGER                                 :: cross
! actual cell
INTEGER                                 :: act_cell

DOUBLE PRECISION                        :: di
INTEGER                                 :: nc

act_cell = package(pack_index)%cell_numb
cross = package(pack_index)%next_cross
IF(cross <= 0 ) THEN
 CALL find_dist(pack_index, nc, di)
 cross = package(pack_index)%next_cross
 write(*,*) 'next_cell_down: act_cell = ', act_cell, ' cross = ', cross
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
   CALL abort()
  END IF
  EXIT
 END IF
END DO

END SUBROUTINE next_cell_down
