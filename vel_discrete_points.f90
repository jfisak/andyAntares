SUBROUTINE vel_discrete_points(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

INTEGER                                 :: testpacket
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
INTEGER                                 :: I, J, K
DOUBLE PRECISION, DIMENSION(7,7)        :: matA
DOUBLE PRECISION, DIMENSION(3, 7)       :: vecB

EXTERNAL                                :: DGETRI
EXTERNAL                                :: DGETRF
EXTERNAL                                :: DGEMN
INTEGER                                 :: info
INTEGER, DIMENSION(size(matA,1))        :: ipiv
DOUBLE PRECISION, DIMENSION(size(matA,1)) :: work

DOUBLE PRECISION                        :: act_coeff
DOUBLE PRECISION, DIMENSION(3, 7)       :: coeffs

testPacket = SIZE(package) - 1
package(testPacket) = package(pack_index)
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
! we will send a testPacket into all six directions
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
 package(testPacket)%dir = directions(cur_dir, :)
 package(testPacket)%next_cross = crossy(cur_dir)
 CALL next_cell_down(testPacket, down_cell)
 IF(dyngrid == 0) THEN
  next_cell = down_cell
 ELSE 
  CALL next_cell_up(testPacket, distances(cur_dir), down_cell, next_cell)
 END IF
 neighbors(cur_dir) = next_cell
END DO

! spherically symmetric model
IF(model_type == 1) THEN
 ! velocity of the current cell
 cur_center = corner + width / 2.0
 cur_vel = cur_vel_norm * cur_center / norm2(cur_center)
 velocity_field(1,:) = cur_vel
 matA(1:2,1) = (/ cur_center(1)**2.0, cur_center(1) /)
 matA(3:4,1) = (/ cur_center(2)**2.0, cur_center(2) /)
 matA(5:7,1) = (/ cur_center(3)**2.0, cur_center(3), 1.D0 /)
 vecB(:, 1) = cur_vel(:)
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
   ! write(*,*) 'vel_discrete_points: n = ', cur_center / norm2(cur_center), &
   !  ' ||v|| = ', cur_vel_norm
   velocity_field(I + 1, :) = cur_vel
  ELSE
   velocity_field(I + 1, :) = velocity_field(1, :)
  END IF
  matA(1:2, I + 1) = (/ cur_center(1)**2.0, cur_center(1) /)
  matA(3:4, I + 1) = (/ cur_center(2)**2.0, cur_center(2) /)
  matA(5:7, I + 1) = (/ cur_center(3)**2.0, cur_center(3), 1.D0 /)
  vecB(:, I + 1) = cur_vel(:)
 END DO

! test for the inverse matrix calculation
! matA(:,1) = (/1, 0, 0, 0, 0, 0, 0 /)
! matA(:,2) = (/0, 1, 0, 0, 0, 0, 0 /)
! matA(:,3) = (/0, 0, 1, 0, 0, 0, 0 /)
! matA(:,4) = (/0, 0, 0, 1, 0, 0, 0 /)
! matA(:,5) = (/0, 0, 0, 0, 1, 0, 0 /)
! matA(:,6) = (/0, 0, 0, 0, 0, 1, 0 /)
! matA(:,7) = (/1, 0, 0, 0, 0, 0, 1 /)

! write(*,*) 'vel_discrete_points: matA = ', matA

CALL DGETRF(7, 7, matA, 7, ipiv, info)
CALL DGETRI(7, matA, 7, ipiv, work, 7, info)

! calculation of coefficients (a ... g)
DO K = 1, 3
 DO I = 1, 7
  act_coeff = 0.0
  DO J = 1, 7
   act_coeff = act_coeff + matA(J, I) * vecB(K, J)
  END DO
  coeffs(K, I) = act_coeff
 END DO
END DO
! write(*,*) 'vel_discrete_points: vecB = ', vecB
DO I = 1, 3
 vel_vec(I) = coeffs(I, 1) * act_pos(1)**2.0 + coeffs(I, 2) * act_pos(1) + &
  coeffs(I, 3) * act_pos(2) + coeffs(I, 4) * act_pos(2) + &
  coeffs(I, 5) * act_pos(3) + coeffs(I, 6) * act_pos(3) + coeffs(I, 7)
END DO
 
write(*,*) 'vel_discrete_points: ||v||/c = ', norm2(vel_vec)/light_speed
! STOP 'vel_discrete_points: testing'

END IF



END SUBROUTINE vel_discrete_points
