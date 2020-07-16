SUBROUTINE vel_discrete_points(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

INTEGER                                 :: dummypackage
INTEGER                                 :: cur_dir
INTEGER, DIMENSION(6,3)                 :: directions
DOUBLE PRECISION, DIMENSION(6)          :: distances
DOUBLE PRECISION, DIMENSION(3)          :: corner, width, act_pos
INTEGER                                 :: next_cell, act_cell, down_cell
INTEGER, DIMENSION(6)                   :: crossy, neighbors
INTEGER                                 :: get_package_model_index, cur_mgi
DOUBLE PRECISION, DIMENSION(3)          :: cur_center, cur_vel, cur_corner, cur_width
DOUBLE PRECISION                        :: cur_vel_norm, cur_vel_norm0
DOUBLE PRECISION, DIMENSION(7,3)        :: velocity_field
INTEGER                                 :: cur_neighbor, cur_nmgi
INTEGER                                 :: I
DOUBLE PRECISION, DIMENSION(7,7)        :: matA

EXTERNAL                                :: DGETRI
EXTERNAL                                :: DGETRF
INTEGER                                 :: info
INTEGER, DIMENSION(size(matA,1))        :: ipiv
DOUBLE PRECISION, DIMENSION(size(matA,1)) :: work

dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
cur_mgi = get_package_model_index(pack_index)
corner = dyn_cell(act_cell)%corner
width = dyn_cell(act_cell)%width
cur_vel_norm0 = model_grid(cur_mgi)%vel

IF(cur_vel_norm == 0.0) THEN
 vel_vec = (/ 0.0, 0.0, 0.0 /)
 RETURN
END IF

! find the neighboring cells
! we will send a dummypackage into all six directions
! x+-, y+-, z+- and find next cells (six cell approximation)

! six possible  directions, distances to the boundaries and crosses
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
 IF(dyngrid == 0) THEN
  next_cell = down_cell
 ELSE 
  CALL next_cell_up(dummypackage, distances(cur_dir), down_cell, next_cell)
 END IF
 neighbors(cur_dir) = next_cell
END DO

! spherically symmetric model
IF(model_type == 1) THEN
 ! velocity of the current cell
 cur_center = corner + width / 2.0
 cur_vel = cur_vel_norm * cur_center / norm2(cur_center)
 velocity_field(1,:) = cur_vel
 matA(1, 1:2) = (/ cur_center(1)**2.0, cur_center(1) /)
 matA(1, 3:4) = (/ cur_center(2)**2.0, cur_center(2) /)
 matA(1, 5:7) = (/ cur_center(3)**2.0, cur_center(3), 1.D0 /)
 ! neighbor velocities
 DO I = 1, 6
  cur_neighbor = neighbors(I)
  IF(cur_neighbor > 0) THEN
   cur_corner = dyn_cell(cur_neighbor)%corner
   cur_width = dyn_cell(cur_neighbor)%width
   cur_center = cur_corner + cur_width / 2.0
   cur_nmgi = dyn_cell(cur_neighbor)%model_index
   cur_vel_norm = model_grid(cur_nmgi)%vel
   cur_vel = cur_vel_norm * cur_center / norm2(cur_center)
   write(*,*) 'vel_discrete_points: n = ', cur_center / norm2(cur_center), &
    ' ||v|| = ', cur_vel_norm
   velocity_field(I + 1, :) = cur_vel
  ELSE
   velocity_field(I + 1, :) = velocity_field(1, :)
  END IF
  matA(I + 1, 1:2) = (/ cur_center(1)**2.0, cur_center(1) /)
  matA(I + 1, 3:4) = (/ cur_center(2)**2.0, cur_center(2) /)
  matA(I + 1, 5:7) = (/ cur_center(3)**2.0, cur_center(3), 1.D0 /)
 END DO

write(*,*) 'vel_discrete_points: matA = ', matA

CALL DGETRF(7, 7, matA, 7, ipiv, info)
CALL DGETRI(7, matA, 7, ipiv, work, 7, info)

write(*,*) 'vel_discrete_points: invmatA = ', matA
 
STOP 'vel_discrete_points: testing'

END IF



END SUBROUTINE vel_discrete_points
