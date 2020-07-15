SUBROUTINE vel_discrete_points(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

INTEGER                                 :: dummypackage
INTEGER                                 :: cur_dir
INTEGER, DIMENSION(6,3)                 :: directions
DOUBLE PRECISION, DIMENSION(6)          :: distances, neighbors
DOUBLE PRECISION, DIMENSION(3)          :: corner, width, act_pos
INTEGER                                 :: next_cross, next_cell, act_cell, down_cell
INTEGER, DIMENSION(6)                   :: crossy 

dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
corner = dyn_cell(act_cell)%corner
width = dyn_cell(act_cell)%width
! find the neighboring cells
! we will send a dummypackage into all six directions
! x+-, y+-, z+- and find next cells (six cell approximation)

! six possible  directions
directions(1,:) = (/ 1, 0, 0 /)
distances(1) = corner(1) + width(1) - act_pos(1)
crossy(1) = posx
directions(2,:) = (/-1, 0, 0 /)
distances(2) = act_pos(1) - corner(1)
crossy(2) = negx
directions(3,:) = (/ 0, 1, 0 /)
distances(3) = corner(2) + width(2) - act_pos(2)
crossy(3) = posy
directions(4,:) = (/ 0,-1, 0 /)
distances(4) = act_pos(2) - corner(2)
crossy(4) = negy
directions(5,:) = (/ 0, 0, 1 /)
distances(5) = corner(3) + width(3) - act_pos(3)
crossy(5) = posz
directions(6,:) = (/ 0, 0,-1 /)
distances(6) = act_pos(3) - corner(3)
crossy(6) = negz

DO cur_dir = 1,6
 package(dummypackage)%dir = directions(cur_dir, :)
 package(dummypackage)%next_cross = crossy(cur_dir)
 CALL next_cell_down(dummypackage, down_cell)
 CALL next_cell_up(dummypackage, distances(cur_dir), down_cell, next_cell)
 neighbors(cur_dir) = next_cell
END DO

STOP 'testing vel_discrete_points'

! spherically symmetric model
IF(model_type == 1) THEN
 
 

END IF



END SUBROUTINE vel_discrete_points
