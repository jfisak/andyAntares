! calculates a factor in equation
! [\mu^2 dv/dr + (1-\mu^2)v/r]
! for several types of velocity fields, concretely
!
! 0 -- homologous approximation
! 1 -- beta law
! 3 -- it is included in the r_kappa_line
! 4 -- calculation via \nabla v
!
! INPUT: pack_index(INT) -- index of packet
!        l_dist(DBLE) -- resonance distance of the next line
!        fr_line(DBLE) -- frequency of the next line
! OUTPUT: roverw(DBLE) -- calculated value of the factor
!
DOUBLE PRECISION FUNCTION roverw(pack_index, l_dist, fr_line)

USE types
USE constants
USE dummypacket
IMPLICIT NONE

DOUBLE PRECISION                               :: R_pos, V_pos, fr_line
DOUBLE PRECISION                                :: l_dist
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: R_pos_vec
DOUBLE PRECISION                                :: costheta
! DOUBLE PRECISION                                :: dV_pos
INTEGER                                         :: dummypack_index, next_cell, pack_index, cur_dummypack

DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: cur_dir
INTEGER, PARAMETER                              :: n_vectors = 3
! DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos_min, pos_pls!, pos_lin
! DOUBLE PRECISION                                :: freq_min, freq_pls
! DOUBLE PRECISION                                :: cmf_min, cmf_pls
! DOUBLE PRECISION                                :: cur_freq_rf
! DOUBLE PRECISION                                :: delta
! DOUBLE PRECISION                                :: deriv
! DOUBLE PRECISION                                :: deriv_min, deriv_pls
! DOUBLE PRECISION, DIMENSION(3)                  :: pos_line
! DOUBLE PRECISION                                :: s_min, s_pls
DOUBLE PRECISION, DIMENSION(n_vectors,const_dimofspace)                :: vel_vectors
LOGICAL                                         :: isposx, isposy, isposz
INTEGER, DIMENSION(n_vectors,const_dimofspace)                         :: directions
INTEGER, DIMENSION(const_dimofspace)                           :: crossy, neighb_cells, mgi_index
DOUBLE PRECISION, DIMENSION(n_vectors)                :: distances
TYPE(dummyphoton)                               :: testPacket

DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: cur_corner, cur_width, act_pos, cur_centre
INTEGER                                         :: down_cell
INTEGER                                         :: cur_pg, cur_xyz_dir

INTEGER                                         :: cur_index_i, cur_index_j, cur_mgi
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: cur_vel, cur_dist
DOUBLE PRECISION                                :: cur_sum, nvn

DOUBLE PRECISION, PARAMETER                     :: large_number = 1.D90
DOUBLE PRECISION                                :: part_1, part_2, part_3
DOUBLE PRECISION, PARAMETER                     :: param_k = 0.4
DOUBLE PRECISION                                :: a_index, b_index, V_star

IF(debug == 100) THEN
 write(*,*) 'roverw: l_dist = ', l_dist
 write(*,*) 'roverw: fr_line = ', fr_line
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cur_dummypack = find_free_index()
dummypack_index = cur_dummypack + SIZE(package)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! write(*,*) 'roverw: calling roverw'

IF(velapprox == 0) THEN
! ONLY FOR TESTING !!!
 IF(sobolev_approximation == 1) THEN
  ROverW = R_inf / V_inf
 END IF ! sobolev_approximation
ELSE IF(velapprox == 2) THEN

 ! according to (10) in Abbot & Lucy (1985)
 ! r
 ! R_pos = norm2(package(dummypack_index)%pos)
 ! ! ||v||
 ! V_pos = V_inf * (1.0 - R_star / R_pos ) ** beta
 ! ! v = (v_x, v_y, v_z)
 ! V_pos_vec = V_pos * package(dummypack_index)%pos / norm2(package(dummypack_index)%pos)
 ! costheta = dot_product(package(dummypack_index)%dir, V_pos_vec) / norm2(V_pos_vec)
 ! ROverW = (V_inf - V_0) / (R_inf - R_star) + 1 / R_pos * (1 - costheta**2.0) * &
 !  (V_0 * R_inf - V_inf * R_star) / (R_inf - R_star)
ELSE IF(velapprox == 1) THEN
 ! according to (10) in Abbot & Lucy (1985)
 ! r
 R_pos = norm2(package(pack_index)%pos)
 ! ||v||
 V_pos = V_inf * (1.0 - R_star / R_pos ) ** beta
 ! write(*,*) 'roverw: V_pos = ', V_pos, ' R_star = ', R_star/R_star, ' R_pos = ', R_pos/R_star, ' V_inf = ', V_inf
 ! v = (v_x, v_y, v_z)
 R_pos_vec = package(pack_index)%pos / norm2(package(pack_index)%pos)
 ! \mu
 costheta = dot_product(package(pack_index)%dir, R_pos_vec)
 ! write(*,*) 'ROverW: dir = ', package(pack_index)%dir, ' V_pos_vec = ', V_pos_vec, ' V_pos = ', V_pos
 ! dv/dr
 ! dV_pos = beta * R_star * V_inf / R_pos**2 * (1.0 - R_star / R_pos)**(beta - 1)
 IF(R_pos <= R_star .or. R_pos > R_inf) THEN
  ROverW = 0.D0
 ELSE
  ROverW = R_pos/V_pos * 1.D0/(costheta**2*((R_star*beta)/(R_pos-R_star)-1.0)+1.0)
 END IF
ELSE IF(velApprox == 5) THEN
 R_pos = norm2(package(pack_index)%pos)
 R_pos_vec = package(pack_index)%pos / norm2(package(pack_index)%pos)
 costheta = dot_product(package(pack_index)%dir, R_pos_vec)
 part_1 = costheta**2*8.0* R_pos/(R_inf - R_star) * cos(8*R_pos/(R_inf - R_star))
 part_2 = (1-costheta**2) * sin(8*R_pos/(R_inf - R_star))
 part_3 = 1/(2 * param_k)
 IF(R_pos <= R_star .or. R_pos > R_inf) THEN
  ROverW = 0.D0
 ELSE
  ROverW = abs(R_pos/(param_k * V_inf) * 1/(part_1 + part_2 + part_3))
 END IF 
!_____________________________________________________________________________________________
! 3D velocity approximation
ELSE IF(velApprox == 3) THEN
 write(*,*) 'roverw: this calculation is implemented into the sbr r_kappa_line'
 STOP 'roverw: exiting'
ELSE IF(velApprox == 10) THEN
 R_pos = norm2(package(pack_index)%pos)
 IF(R_pos <= R_star .or. R_pos > R_inf) THEN
  ROverW = 0.D0
 ELSE
  R_pos = norm2(package(pack_index)%pos)
  R_pos_vec = package(pack_index)%pos / norm2(package(pack_index)%pos)
  costheta = dot_product(package(pack_index)%dir, R_pos_vec)
  V_star = R_star/R_inf * V_inf
  a_index = - (V_inf - V_star)/(R_inf - R_star)
  b_index = (V_inf * R_inf - V_star * R_star)/(R_inf - R_star)
  ROverW = (abs(a_index + (1-costheta**2)*b_index/R_pos))**(-1)
 END IF 
 

! calculating nabla v during the packet propagation
ELSE IF(velApprox == 4) THEN

 act_pos = package(pack_index)%pos
 cur_pg = package(pack_index)%cell_numb
 cur_dir = package(pack_index)%dir
 cur_corner = dyn_cell(cur_pg)%corner
 cur_width = dyn_cell(cur_pg)%width
 cur_centre = cur_corner + cur_width/2.0
 IF(vel_propgrid) THEN
  cur_vel = dyn_cell(cur_pg)%vec_vel
 ELSE IF(vel_modgrid) THEN
  cur_mgi = dyn_cell(cur_pg)%model_index
  cur_vel = model_grid(cur_mgi)%vec_vel
 END IF

 ! a calculation of a velocity gradient
 ! velocities in three different directions
 ! according to this quantities we will choose the neighbouring cells
 ! x
 IF(package(pack_index)%dir(ind_x) > 0) THEN
  isposx = .true.
 ELSE
  isposx = .false.
 END IF

 ! y
 IF(package(pack_index)%dir(ind_y) > 0) THEN
  isposy = .true.
 ELSE
  isposy = .false.
 END IF

 ! z
 IF(package(pack_index)%dir(ind_z) > 0) THEN
  isposz = .true.
 ELSE
  isposz = .false.
 END IF
 
 ! find the neighboring cells
 ! we will send a testPacket into all six directions
 ! x+-, y+-, z+- and find next cells (six cell approximation)
 
 ! six possible  directions, distances to the boundaries and crosses
 IF(isposx) THEN
  directions(ind_x,:) = (/ 1, 0, 0 /)
  distances(ind_x) = cur_corner(ind_x) + cur_width(ind_x) - act_pos(ind_x)
  crossy(ind_x) = posx
 ELSE
  directions(ind_x,:) = (/-1, 0, 0 /)
  distances(ind_x) = act_pos(ind_x) - cur_corner(ind_x)
  crossy(ind_x) = negx
 END IF

 IF(isposy) THEN
  directions(ind_y,:) = (/ 0, 1, 0 /)
  distances(ind_y) = cur_corner(ind_y) + cur_width(ind_y) - act_pos(ind_y)
  crossy(ind_y) = posy
 ELSE
  directions(ind_y,:) = (/ 0,-1, 0 /)
  distances(ind_y) = act_pos(ind_y) - cur_corner(ind_y)
  crossy(ind_y) = negy
 END IF

 IF(isposz) THEN
  directions(ind_z,:) = (/ 0, 0,  1 /)
  distances(ind_z) = cur_corner(ind_z) + cur_width(ind_z) - act_pos(ind_z)
  crossy(ind_z) = posz
 ELSE
  directions(ind_z,:) = (/ 0, 0, -1 /)
  distances(ind_z) = act_pos(ind_z) - cur_corner(ind_z)
  crossy(ind_z) = negz
 END IF

 DO cur_xyz_dir = 1, const_dimofspace
  testPacket%dir = directions(cur_xyz_dir, :)
  testPacket%next_cross = crossy(cur_xyz_dir)
  testPacket%cell_numb = package(pack_index)%cell_numb
  testPacket%pos = package(pack_index)%pos
  ! write(*,*) 'roverw: dummypack_index = ', dummypack_index
  dummypackage(cur_dummypack) = testPacket
  CALL next_cell_down(dummypack_index, down_cell)
  IF(dyngrid == 0) THEN
   next_cell = down_cell
  ELSE 
   CALL next_cell_up(dummypack_index, distances(cur_xyz_dir), down_cell, next_cell)
  END IF
  ! write(*,*) 'roverw: next_cell = ', next_cell
  
  neighb_cells(cur_xyz_dir) = next_cell
  cur_dist(cur_xyz_dir) = (dyn_cell(cur_pg)%width(cur_xyz_dir) + dyn_cell(next_cell)%width(cur_xyz_dir))/2.0
 END DO

write(*,*) 'roverw: vel_modgrid = ', vel_modgrid, ' vel_propgrid = ', vel_propgrid
IF(velApprox == 3 .and. vel_modgrid) THEN
 IF(model_type == 3) THEN
  mgi_index(ind_x) = dyn_cell(neighb_cells(ind_x))%model_index
  vel_vectors(ind_x,:) = model_grid(mgi_index(ind_x))%vec_vel
  mgi_index(ind_y) = dyn_cell(neighb_cells(ind_y))%model_index
  vel_vectors(ind_y,:) = model_grid(mgi_index(ind_y))%vec_vel
  mgi_index(ind_z) = dyn_cell(neighb_cells(ind_z))%model_index
  vel_vectors(ind_z,:) = model_grid(mgi_index(ind_z))%vec_vel
 ELSE IF(model_type == 3) THEN
  STOP 'vel_discrete_points: this type of intput is not supported yet'
 END IF
ELSE IF(velApprox == 4 .or. vel_propgrid) THEN ! velocity is pre-calculated in each propagation cell
 vel_vectors(ind_x,:) = dyn_cell(neighb_cells(ind_x))%vec_vel
 vel_vectors(ind_y,:) = dyn_cell(neighb_cells(ind_y))%vec_vel
 vel_vectors(ind_z,:) = dyn_cell(neighb_cells(ind_z))%vec_vel
END IF

! the calculation of the velocity gradient
DO cur_index_i = 1, const_dimofspace ! index of a derivative
 DO cur_index_j = 1, const_dimofspace ! index of a vector v
  vel_vectors(cur_index_i,cur_index_j) = (cur_vel(cur_index_j) - vel_vectors(cur_index_i, cur_index_j))/cur_dist(cur_index_i)
  ! write(*,*) 'roverw: vel_0 = ', cur_vel(cur_index_j), ' vel_j = ', vel_vectors(cur_index_i, cur_index_j)
  ! write(*,*) 'roverw: dist = ', cur_dist(cur_index_i)
  ! write(*,*) 'roverw: vel_vectors( ', cur_index_i, ', ', cur_index_j, ' ) = ', vel_vectors(cur_index_i,cur_index_j)
 END DO
END DO

nvn = 0.D0
! a multiplication with direction vectors
DO cur_index_i = 1, const_dimofspace
 cur_sum = 0.D0
 DO cur_index_j = 1, const_dimofspace
  cur_sum = cur_sum + cur_dir(cur_index_j) * vel_vectors(cur_index_i, cur_index_j)
 END DO
 nvn = nvn + cur_sum * cur_dir(cur_index_i)
END DO

! and finally the absolute value
nvn = abs(nvn)
! write(*,*) 'roverw: nvn = ', nvn

IF(nvn == 0.D0) THEN
 roverw = large_number
ELSE
 roverw = 1/nvn
END IF
ELSE
 write(*,*) 'roverw: velapprox = ', velapprox, ' is not a valid choice'
 STOP
END IF

CALL deactivate_dummy_packet(cur_dummypack)


END FUNCTION
