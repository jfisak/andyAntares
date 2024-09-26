! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).
!
! INPUT: pack_index, INT -- index of a packet
! OUTPUT: dist, DBLE -- distance to the closest boundary of the propGrid cell
!         next_cell, INT -- index of the next propGrid cell
!
SUBROUTINE boundary3(pack_index, dist, next_cell) 

USE types
USE constants

IMPLICIT NONE    

INTEGER                         :: pack_index, next_cell, n_cell
DOUBLE PRECISION                :: dist
! dynamic cell variables
INTEGER                         :: act_cell, next_cross

! calculation of a distance from the basic cell
! firstly we have to know which basic cell photon occupies
!   Number of the current cell
act_cell = package(pack_index)%cell_numb
! now we are computing the nearest distance to the actuall dynamic cell
CALL bound_dist(pack_index, act_cell, dist)
next_cross = package(pack_index)%next_cross
IF(next_cross <= 6) THEN
 ! now we look for the next cell given by indexes
 CALL next_cell_down(pack_index, n_cell)
 ! position of the point
 !cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist
 if(dyngrid /= 0) CALL next_cell_up(pack_index, dist, n_cell, next_cell)
 if(dyngrid == 0) next_cell = n_cell
ELSE
 next_cell = act_cell
END IF
! write(*,*) 'boundary3: next_cell = ', next_cell
IF(next_cell > SIZE(dyn_cell)) THEN
 write(*,*) 'boundary3: pack_index = ', pack_index, ' dist = ', dist, &
  ' next_cell = ', next_cell, ' is larger than the size of dyn_cell'
  CALL abort()
END IF

END SUBROUTINE boundary3
