SUBROUTINE boundary3(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

USE types

IMPLICIT NONE    

INTEGER                         :: pack_index, next_cell, n_cell
DOUBLE PRECISION                :: dist
! dynamic cell variables
INTEGER                         :: act_cell

! calculation of a distance from the basic cell
! firstly we have to know which basic cell photon occupies
!   Number of the current cell
act_cell = package(pack_index)%cell_numb
!  print*, 'boundary3: pack index = ', pack_index, ' actCell = ', actCell
!  print*, 'boundary3: basic_cell_numb = ', basic_cell_numb
! now we are computing the nearest distance to the actuall dynamic cell
CALL find_dist(pack_index, act_cell, dist)
! now we look for the next cell given by indexes
CALL next_cell_down(pack_index, n_cell)
! position of the point
!cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist
if(dyngrid /= 0) CALL next_cell_up(pack_index, dist, n_cell, next_cell)
if(dyngrid == 0) next_cell = n_cell
! write(*,*) 'boundary3: pack_index = ', pack_index, ' dist = ', dist, &
!  ' next_cell = ', next_cell
IF(next_cell > SIZE(dyn_cell)) THEN
 write(*,*) 'boundary3: pack_index = ', pack_index, ' dist = ', dist, &
  ' next_cell = ', next_cell, ' is larger than the size of dyn_cell'
  CALL abort()
END IF

END SUBROUTINE boundary3
