! this sbr interpolates a velocity in the defined position in the propGrid
! it uses the trilinear interpolation
!
! INPUT: pack_index -- an index of a packet
!
! OUTPUT: vel_vec -- calculated velocity vector
!
! RETURN point: 1x
!
SUBROUTINE vel_discrete_points_interpolation(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: vel_vec

TYPE(photon)                            :: testPacket
INTEGER                                 :: dummyPacket
INTEGER                                 :: cur_dir
INTEGER, DIMENSION(6,const_dimofspace)                 :: directions
DOUBLE PRECISION, DIMENSION(6)          :: distances
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: corner, width, act_pos
INTEGER                                 :: next_cell, act_cell, down_cell
INTEGER, DIMENSION(6)                   :: crossy, neighbors
INTEGER                                 :: get_package_model_index, cur_mgi
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cur_center, cur_vel, cur_corner, cur_width
DOUBLE PRECISION                        :: cur_vel_norm, cur_vel_norm0
DOUBLE PRECISION, DIMENSION(7,const_dimofspace)        :: velocity_field
INTEGER                                 :: cur_neighbor, cur_nmgi
INTEGER                                 :: ind_I
DOUBLE PRECISION, DIMENSION(7,7)        :: matA
DOUBLE PRECISION, DIMENSION(const_dimofspace, 7)       :: vecB

EXTERNAL                                :: DGETRI
EXTERNAL                                :: DGETRF
EXTERNAL                                :: DGEMN
INTEGER                                 :: info
INTEGER, DIMENSION(size(matA,1))        :: ipiv
DOUBLE PRECISION, DIMENSION(size(matA,1)) :: work

DOUBLE PRECISION                        :: act_coeff
DOUBLE PRECISION, DIMENSION(const_dimofspace, 7)       :: coeffs

DOUBLE PRECISION                        :: sumx, sumxsq, sumy, sumysq, sumxy
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: vel_m, vel_p, vel_0
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: pos_m, pos_p, pos_0
INTEGER                                 :: cur_dim, n_index
INTEGER                                 :: cur_neighbor_p, cur_neighbor_m
INTEGER                                 :: n_points

INTEGER                                 :: mgi_p, mgi_m
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: inda, indb
DOUBLE PRECISION                        :: vel_ana

testPacket = package(pack_index)
dummyPacket = SIZE(package)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
cur_mgi = get_package_model_index(pack_index)
corner = dyn_cell(act_cell)%corner
upcorner = dyn_cell(act_cell)%upcorner
cur_center = corner + width / 2.0
width = dyn_cell(act_cell)%width
cur_vel_norm0 = model_grid(cur_mgi)%vel

IF(cur_vel_norm == 0.0) THEN
 vel_vec = (/ 0.0, 0.0, 0.0 /)
 ! RETURN point
 RETURN
END IF

! find the neighboring cells
! we will send a testPacket into all six directions
! x+-, y+-, z+- and find next cells (six cell approximation)

! six possible  directions, distances to the boundaries and crosses
directions(posx,:) = (/ 1, 0, 0 /)
distances(posx) = upcorner(ind_x) - act_pos(ind_x)
crossy(posx) = posx
directions(negx,:) = (/-1, 0, 0 /)
distances(negx) = act_pos(ind_x) - corner(ind_x)
crossy(negx) = negx
directions(posy,:) = (/ 0, 1, 0 /)
distances(posy) = upcorner(ind_y)- act_pos(ind_y)
crossy(posy) = posy
directions(negy,:) = (/ 0,-1, 0 /)
distances(negy) = act_pos(ind_y) - corner(ind_y)
crossy(negy) = negy
directions(posz,:) = (/ 0, 0, 1 /)
distances(posz) = upcorner(ind_z) - act_pos(ind_z)
crossy(posz) = posz
directions(negz,:) = (/ 0, 0,-1 /)
distances(negz) = act_pos(ind_z) - corner(ind_z)
crossy(negz) = negz

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
 IF(velApprox == 3 .and. vel_modgrid) THEN
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
  ELSE IF(model_type == 3) THEN
   STOP 'vel_discrete_points: this type of intput is not supported yet'
  END IF
 ELSE IF(velApprox == 4 .or. vel_propgrid) THEN ! velocity is pre-calculated in each propagation cell
  vel_m = dyn_cell(cur_neighbor_m)%vec_vel
  vel_p = dyn_cell(cur_neighbor_p)%vec_vel
  vel_0 = dyn_cell(act_cell)%vec_vel
 END IF
 ! sums
 sumx = pos_m(cur_dim) + pos_0(cur_dim) + pos_p(cur_dim)
 sumxsq = pos_m(cur_dim)**2 + pos_0(cur_dim)**2 + pos_p(cur_dim)**2
 sumy = vel_m(cur_dim) + vel_0(cur_dim) + vel_p(cur_dim)
 sumysq = vel_m(cur_dim)**2 + vel_0(cur_dim)**2 + vel_p(cur_dim)**2
 sumxy = pos_m(cur_dim) * vel_m(cur_dim) + pos_0(cur_dim) * vel_0(cur_dim) + &
  & pos_p(cur_dim) * vel_p(cur_dim)
 write(*,*) 'vel_discrete_points: sumx = ', sumx, ' sumy = ', sumy

 inda(cur_dim) = (n_points * sumxy - sumx * sumy) / (n_points * sumxsq - sumx**2)
 indb(cur_dim) = (sumxsq * sumy - sumx * sumxy)   / (n_points * sumxsq - sumx**2)
 n_index = n_index + 2
END DO

DO ind_I = 1, 3
 vel_vec(ind_I) = inda(ind_I) * act_pos(ind_I) + indb(ind_I)
END DO

vel_ana = V_inf/R_inf * norm2(package(pack_index)%pos)
write(*,*) 'vel_discrete_points: ||v||/v_an = ', norm2(vel_vec)/vel_ana




END SUBROUTINE vel_discrete_points_interpolation
