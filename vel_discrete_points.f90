SUBROUTINE vel_discrete_points(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

TYPE(photon)                            :: testPacket
INTEGER                                 :: dummyPacket
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

DOUBLE PRECISION                        :: sumx, sumxsq, sumy, sumysq, sumxy
DOUBLE PRECISION, DIMENSION(3)          :: vel_m, vel_p, vel_0
DOUBLE PRECISION, DIMENSION(3)          :: pos_m, pos_p, pos_0
INTEGER                                 :: cur_dim, n_index
INTEGER                                 :: cur_neighbor_p, cur_neighbor_m
INTEGER                                 :: n_points

INTEGER                                 :: mgi_p, mgi_m
DOUBLE PRECISION, DIMENSION(3)          :: inda, indb
DOUBLE PRECISION                        :: vel_ana

testPacket = package(pack_index)
dummyPacket = SIZE(package)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
cur_mgi = get_package_model_index(pack_index)
corner = dyn_cell(act_cell)%corner
cur_center = corner + width / 2.0
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

! searching the neighbor cells
DO cur_dir = 1,6
 testPacket%dir = directions(cur_dir, :)
 testPacket%next_cross = crossy(cur_dir)
 package(dummyPacket) = testPacket
 CALL next_cell_down(dummyPacket, down_cell)
 IF(dyngrid == 0) THEN
  next_cell = down_cell
 ELSE 
  CALL next_cell_up(dummyPacket, distances(cur_dir), down_cell, next_cell)
 END IF
 neighbors(cur_dir) = next_cell
END DO

! spherically symmetric model
! IF(model_type == 1) THEN
!  ! velocity of the current cell
!  cur_vel = cur_vel_norm * cur_center / norm2(cur_center)
!  velocity_field(1,:) = cur_vel
! END IF

n_points = 3
n_index = 1
DO cur_dim = 1, 3
 cur_neighbor_p = neighbors(n_index)
 cur_neighbor_m = neighbors(n_index + 1)
 ! positions

 IF (cur_neighbor_m > 0) THEN
  pos_m = dyn_cell(cur_neighbor_m)%corner + dyn_cell(cur_neighbor_m)%width/2.0
  mgi_m = dyn_cell(cur_neighbor_m)%model_index
 ! if it is out of the grid we will set it to the actuall cell boundary
 ELSE 
  pos_m = act_pos
  pos_m(cur_dim) = corner(cur_dim)
  mgi_m = cur_mgi
 END IF
 IF (cur_neighbor_p > 0) THEN
  pos_p = dyn_cell(cur_neighbor_p)%corner + dyn_cell(cur_neighbor_p)%width/2.0
  mgi_p = dyn_cell(cur_neighbor_p)%model_index
 ELSE 
  pos_m = act_pos
  pos_m(cur_dim) = corner(cur_dim) + width(cur_dim)
  mgi_p = cur_mgi
 END IF
 pos_0 = cur_center
 ! velocities
 IF(model_type == 1) THEN
  vel_m = model_grid(mgi_m)%vel * pos_m / norm2(pos_m)
  vel_p = model_grid(mgi_p)%vel * pos_p / norm2(pos_p)
  vel_0 = model_grid(cur_mgi)%vel * pos_0 / norm2(pos_0)
  write(*,*) 'vel_discrete_points: vel_0 = ', vel_0, ' pos_0 = ', pos_0
 ELSE IF(model_type == 2) THEN
  mgi_m = dyn_cell(cur_neighbor_m)%model_index
  vel_m = model_grid(mgi_m)%velocity
  mgi_p = dyn_cell(cur_neighbor_p)%model_index
  vel_p = model_grid(mgi_p)%velocity
  vel_0 = model_grid(cur_mgi)%velocity
 END IF
 ! sums
 sumx = pos_m(cur_dim) + pos_0(cur_dim) + pos_p(cur_dim)
 sumxsq = pos_m(cur_dim)**2 + pos_0(cur_dim)**2 + pos_p(cur_dim)**2
 sumy = vel_m(cur_dim) + vel_0(cur_dim) + vel_p(cur_dim)
 sumysq = vel_m(cur_dim)**2 + vel_0(cur_dim)**2 + vel_p(cur_dim)**2
 sumxy = pos_m(cur_dim) * vel_m(cur_dim) + pos_0(cur_dim) * vel_0(cur_dim) + &
  pos_p(cur_dim) * vel_p(cur_dim)
 write(*,*) 'vel_discrete_points: sumx = ', sumx, ' sumy = ', sumy

 inda(cur_dim) = (n_points * sumxy - sumx * sumy)/(n_points * sumxsq - sumx**2)
 indb(cur_dim) = (sumxsq * sumy - sumx * sumxy)/(n_points * sumxsq - sumx**2)
 n_index = n_index + 2
END DO

DO I = 1, 3
 vel_vec(I) = inda(I) * act_pos(I) + indb(I)
END DO

vel_ana = V_inf/R_inf * norm2(package(pack_index)%pos)
write(*,*) 'vel_discrete_points: ||v||/v_an = ', norm2(vel_vec)/vel_ana
!  matA(1:2,1) = (/ cur_center(1)**2.0, cur_center(1) /)
!  matA(3:4,1) = (/ cur_center(2)**2.0, cur_center(2) /)
!  matA(5:7,1) = (/ cur_center(3)**2.0, cur_center(3), 1.D0 /)
!  vecB(:, 1) = cur_vel(:)
!  ! neighbor velocities
!  DO I = 1, 6
!   cur_neighbor = neighbors(I)
!   IF(cur_neighbor > 0) THEN
!    cur_corner = dyn_cell(cur_neighbor)%corner
!    cur_width = dyn_cell(cur_neighbor)%width
!    cur_center = cur_corner + cur_width / 2.0
!    cur_nmgi = dyn_cell(cur_neighbor)%model_index
!    cur_vel_norm = model_grid(cur_nmgi)%vel
!    cur_vel = cur_vel_norm * cur_center / norm2(cur_center)
!    ! write(*,*) 'vel_discrete_points: n = ', cur_center / norm2(cur_center), &
!    !  ' ||v|| = ', cur_vel_norm
!    velocity_field(I + 1, :) = cur_vel
!   ELSE
!    velocity_field(I + 1, :) = velocity_field(1, :)
!   END IF
!   matA(1:2, I + 1) = (/ cur_center(1)**2.0, cur_center(1) /)
!   matA(3:4, I + 1) = (/ cur_center(2)**2.0, cur_center(2) /)
!   matA(5:7, I + 1) = (/ cur_center(3)**2.0, cur_center(3), 1.D0 /)
!   vecB(:, I + 1) = cur_vel(:)
!  END DO

! test for the inverse matrix calculation
! matA(:,1) = (/1, 0, 0, 0, 0, 0, 0 /)
! matA(:,2) = (/0, 1, 0, 0, 0, 0, 0 /)
! matA(:,3) = (/0, 0, 1, 0, 0, 0, 0 /)
! matA(:,4) = (/0, 0, 0, 1, 0, 0, 0 /)
! matA(:,5) = (/0, 0, 0, 0, 1, 0, 0 /)
! matA(:,6) = (/0, 0, 0, 0, 0, 1, 0 /)
! matA(:,7) = (/1, 0, 0, 0, 0, 0, 1 /)

! write(*,*) 'vel_discrete_points: matA = ', matA

! CALL DGETRF(7, 7, matA, 7, ipiv, info)
! CALL DGETRI(7, matA, 7, ipiv, work, 7, info)
! 
! ! calculation of coefficients (a ... g)
! DO K = 1, 3
!  DO I = 1, 7
!   act_coeff = 0.0
!   DO J = 1, 7
!    act_coeff = act_coeff + matA(J, I) * vecB(K, J)
!   END DO
!   coeffs(K, I) = act_coeff
!  END DO
! END DO
! ! write(*,*) 'vel_discrete_points: vecB = ', vecB
! DO I = 1, 3
!  vel_vec(I) = coeffs(I, 1) * act_pos(1)**2.0 + coeffs(I, 2) * act_pos(1) + &
!   coeffs(I, 3) * act_pos(2) + coeffs(I, 4) * act_pos(2) + &
!   coeffs(I, 5) * act_pos(3) + coeffs(I, 6) * act_pos(3) + coeffs(I, 7)
! END DO
 
! STOP 'vel_discrete_points: testing'




END SUBROUTINE vel_discrete_points
