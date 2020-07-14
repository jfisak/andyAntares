SUBROUTINE vel_discrete_points(pack_index, vel_vec)

IMPLICIT NONE

USE types

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

INTEGER                                 :: dummypackage
INTEGER, DIMENSION(6,3)                 :: directions
DOUBLE PRECISION, DIMENSION(6)          :: distances

dummypackage = SIZE(package)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
corner = dyn_cell(act_cell)%corner
width = dyn_cell(act_cell)%width
! find the neighboring cells
! we will send a dummypackage into all six directions
! x+-, y+-, z+- and find next cells (six cell approximation)
package(dummypackage)

! six possible  directions
directions(1,:) = (/ 1, 0, 0 /)
distances(1) = 
directions(2,:) = (/-1, 0, 0 /)
directions(3,:) = (/ 0, 1, 0 /)
directions(4,:) = (/ 0,-1, 0 /)
directions(5,:) = (/ 0, 0, 1 /)
directions(6,:) = (/ 0, 0,-1 /)

DO cur_dir = 1,6
 package(dummypackage)%dir = directions(cur_dir, :)
 
END DO

! spherically symmetric model
IF(model_type == 1) THEN
 
 

END IF



END SUBROUTINE vel_discrete_points
