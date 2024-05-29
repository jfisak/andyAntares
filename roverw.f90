DOUBLE PRECISION FUNCTION roverw(pack_index, l_dist, fr_line)

USE types
USE constants
USE dummypacket
IMPLICIT NONE

DOUBLE PRECISION                               :: R_pos, V_pos, fr_line
DOUBLE PRECISION                                :: l_dist
DOUBLE PRECISION, DIMENSION(3)                  :: V_pos_vec
DOUBLE PRECISION                                :: costheta
DOUBLE PRECISION                                :: dV_pos
DOUBLE PRECISION                                :: cell_dist
INTEGER                                         :: dummypack_index, dummypack_index0, next_cell, pack_index

DOUBLE PRECISION, DIMENSION(3)                  :: cur_pos, cur_dir
DOUBLE PRECISION, DIMENSION(3)                  :: pos_min, pos_pls!, pos_lin
! DOUBLE PRECISION                                :: freq_min, freq_pls
DOUBLE PRECISION                                :: cmf_min, cmf_pls
DOUBLE PRECISION                                :: cur_freq_rf
! DOUBLE PRECISION, PARAMETER                     :: delta=1.E0
! DOUBLE PRECISION                                :: delta
DOUBLE PRECISION                                :: deriv
! DOUBLE PRECISION                                :: deriv_min, deriv_pls
! DOUBLE PRECISION, DIMENSION(3)                  :: pos_line
DOUBLE PRECISION                                :: s_min, s_pls
DOUBLE PRECISION, DIMENSION(3,3)                :: deriVel, vel_vectors
LOGICAL                                         :: isposx, isposy, isposz
INTEGER, DIMENSION(3,3)                         :: directions
INTEGER, DIMENSION(3)                           :: crossy, neighb_cells, mgi_index
DOUBLE PRECISION, DIMENSION(3)                :: distances
TYPE(dummyphoton)                               :: testPacket

DOUBLE PRECISION, DIMENSION(3)                  :: cur_corner, cur_width, act_pos, cur_centre
INTEGER                                         :: down_cell
INTEGER                                         :: cur_pg, cur_xyz_dir
DOUBLE PRECISION, DIMENSION(3,3)                :: tensor_divj

INTEGER                                         :: cur_index_i, cur_index_j, cur_mgi
DOUBLE PRECISION, DIMENSION(3)                  :: cur_vel, cur_dist
DOUBLE PRECISION                                :: cur_sum, nvn

DOUBLE PRECISION, PARAMETER                     :: large_number = 1.D90

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! toto budu muset ještě změnit
! DOČASNÉ ŘEŠENÍ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
dummypack_index = 1
dummypack_index0 = SIZE(package) - 3 + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF(velapprox == 0) THEN
! ONLY FOR TESTING !!!
 IF(sobolev_approximation == 1) THEN
  ROverW = R_inf / V_inf
 END IF ! sobolev_approximation
ELSE IF(velapprox == 2) THEN

 ! according to (10) in Abbot & Lucy (1985)
 ! r
 R_pos = norm2(package(dummypack_index)%pos)
 ! ||v||
 V_pos = V_inf * (1.0 - R_star / R_pos ) ** beta
 ! v = (v_x, v_y, v_z)
 V_pos_vec = V_pos * package(dummypack_index)%pos / norm2(package(dummypack_index)%pos)
 costheta = dot_product(package(dummypack_index)%dir, V_pos_vec) / norm2(V_pos_vec)
 ROverW = (V_inf - V_0) / (R_inf - R_star) + 1 / R_pos * (1 - costheta**2.0) * &
  (V_0 * R_inf - V_inf * R_star) / (R_inf - R_star)
ELSE IF(velapprox == 1) THEN
 CALL boundary3(pack_index, cell_dist, next_cell)
 ! according to (10) in Abbot & Lucy (1985)
 ! r
 R_pos = norm2(package(pack_index)%pos)
 ! ||v||
 V_pos = V_inf * (1.0 - R_star / R_pos ) ** beta
 ! v = (v_x, v_y, v_z)
 V_pos_vec = V_pos * package(pack_index)%pos / norm2(package(pack_index)%pos)
 ! \mu
 costheta = dot_product(package(pack_index)%dir, V_pos_vec) / norm2(V_pos_vec)
 ! dv/dr
 dV_pos = beta * R_star * V_inf / R_pos**2 * (1.0 - R_star / R_pos)**(beta-1)
 ROverW = 1.0 / (costheta**2.0 * dV_pos + (1.0 - costheta**2.0)* V_pos / R_pos)
 ! actirrates%Lline(I) = low_pop * Blu * h * light_speed * &
 !  ROverW / (4.0 * pi) * corrFactor 
!_____________________________________________________________________________________________
! 3D velocity approximation
ELSE IF(velapprox == 3) THEN
 ! we will have to find CMF frequencies at two points, the middle location is the Sobolev point
 IF(sobolev_approximation == 1) THEN
  
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
  IF(package(pack_index)%dir(1) > 0) THEN
   isposx = .true.
  ELSE
   isposx = .false.
  END IF

  ! y
  IF(package(pack_index)%dir(2) > 0) THEN
   isposy = .true.
  ELSE
   isposy = .false.
  END IF

  ! z
  IF(package(pack_index)%dir(3) > 0) THEN
   isposz = .true.
  ELSE
   isposz = .false.
  END IF
  
  ! find the neighboring cells
  ! we will send a testPacket into all six directions
  ! x+-, y+-, z+- and find next cells (six cell approximation)
  
  ! six possible  directions, distances to the boundaries and crosses
  IF(isposx) THEN
   directions(1,:) = (/ 1, 0, 0 /)
   distances(1) = cur_corner(1) + cur_width(1) - act_pos(1)
   crossy(1) = posx
  ELSE
   directions(1,:) = (/-1, 0, 0 /)
   distances(1) = act_pos(1) - cur_corner(1)
   crossy(1) = negx
  END IF

  IF(isposy) THEN
   directions(2,:) = (/ 0, 1, 0 /)
   distances(2) = cur_corner(2) + cur_width(2) - act_pos(2)
   crossy(2) = posy
  ELSE
   directions(2,:) = (/ 0,-1, 0 /)
   distances(2) = act_pos(2) - cur_corner(2)
   crossy(2) = negy
  END IF

  IF(isposz) THEN
   directions(3,:) = (/ 0, 0,  1 /)
   distances(3) = cur_corner(3) + cur_width(3) - act_pos(3)
   crossy(3) = posz
  ELSE
   directions(3,:) = (/ 0, 0, -1 /)
   distances(3) = act_pos(3) - cur_corner(3)
   crossy(3) = negz
  END IF

  DO cur_xyz_dir = 1,3
   testPacket%dir = directions(cur_xyz_dir, :)
   testPacket%next_cross = crossy(cur_xyz_dir)
   testPacket%cell_numb = package(pack_index)%cell_numb
   testPacket%pos = package(pack_index)%pos
   ! write(*,*) 'roverw: dummypack_index = ', dummypack_index
   dummypackage(dummypack_index) = testPacket
   CALL next_cell_down(SIZE(package) + 10 + dummypack_index, down_cell)
   IF(dyngrid == 0) THEN
    next_cell = down_cell
   ELSE 
    CALL next_cell_up(SIZE(package) + 10 + dummypack_index, distances(cur_xyz_dir), down_cell, next_cell)
   END IF
   ! write(*,*) 'roverw: next_cell = ', next_cell
   
   neighb_cells(cur_xyz_dir) = next_cell
   cur_dist(cur_xyz_dir) = (dyn_cell(cur_pg)%width(cur_xyz_dir) + dyn_cell(next_cell)%width(cur_xyz_dir))/2.0
  END DO
  
 ! write(*,*) 'roverw: vel_modgrid = ', vel_modgrid, ' vel_propgrid = ', vel_propgrid
 IF(velApprox == 3 .and. vel_modgrid) THEN
  IF(model_type == 3) THEN
   mgi_index(1) = dyn_cell(neighb_cells(1))%model_index
   vel_vectors(1,:) = model_grid(mgi_index(1))%vec_vel
   mgi_index(2) = dyn_cell(neighb_cells(2))%model_index
   vel_vectors(2,:) = model_grid(mgi_index(2))%vec_vel
   mgi_index(3) = dyn_cell(neighb_cells(3))%model_index
   vel_vectors(3,:) = model_grid(mgi_index(3))%vec_vel
  ELSE IF(model_type == 3) THEN
   STOP 'vel_discrete_points: this type of intput is not supported yet'
  END IF
 ELSE IF(velApprox == 4 .or. vel_propgrid) THEN ! velocity is pre-calculated in each propagation cell
  vel_vectors(1,:) = dyn_cell(neighb_cells(1))%vec_vel
  vel_vectors(2,:) = dyn_cell(neighb_cells(2))%vec_vel
  vel_vectors(3,:) = dyn_cell(neighb_cells(3))%vec_vel
 END IF

 ! the calculation of the velocity gradient
 DO cur_index_i = 1,3 ! index of a derivative
  DO cur_index_j = 1,3 ! index of a vector v
   vel_vectors(cur_index_i,cur_index_j) = (cur_vel(cur_index_j) - vel_vectors(cur_index_i, cur_index_j))/cur_dist(cur_index_i)
   ! write(*,*) 'roverw: vel_0 = ', cur_vel(cur_index_j), ' vel_j = ', vel_vectors(cur_index_i, cur_index_j)
   ! write(*,*) 'roverw: dist = ', cur_dist(cur_index_i)
   ! write(*,*) 'roverw: vel_vectors( ', cur_index_i, ', ', cur_index_j, ' ) = ', vel_vectors(cur_index_i,cur_index_j)
  END DO
 END DO

 nvn = 0.D0
 ! a multiplication with direction vectors
 DO cur_index_i = 1,3
  cur_sum = 0.D0
  DO cur_index_j = 1,3
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


 !  cur_pos = package(pack_index)%pos
 !  cur_dir = package(pack_index)%dir
 !  cur_freq_rf = package(pack_index)%freq_rf

 !  s_min = l_dist - delta
 !  s_pls = l_dist + delta

 !  pos_min = cur_pos + cur_dir * s_min
 !  pos_pls = cur_pos + cur_dir * s_pls

 !  ! write(*,*) 'roverw: pos_min = ', norm2(pos_min), ' pos_pls = ', norm2(pos_pls)
 !  ! write(*,*) 'roverw: pos/s_min = ', norm2(cur_pos)/s_min

 !  CALL cmf_freq(pack_index, pos_min, cur_freq_rf, cmf_min)
 !  CALL cmf_freq(pack_index, pos_pls, cur_freq_rf, cmf_pls)

 !  IF(cmf_min == cmf_pls) THEN
 !   write(*,*) 'roverw: l_dist = ', l_dist/R_star
 !   STOP 'roverw: cmf_min == cmf_pls'
 !  END IF
 !  ! write(*,*) 'roverw: p+ - p- = ', pos_pls - pos_min
 !  ! write(*,*) 'roverw: pos_min = ', pos_min, ' pos_pls = ', pos_pls
 !  ! write(*,*) 'roverw: f+ - f- = ', cmf_pls - cmf_min

 !  deriv = (s_pls - s_min)/(cmf_pls - cmf_min)
 !  write(*,*) 'roverw deriv = ', deriv
 !  write(*,*) 'roverw: f-, f0, f+ = ', cmf_min, fr_line, cmf_pls
 !  write(*,*) 'roverw: s-, s0, s+ = ', s_min, s_pls

 !  roverw = deriv
 !  write(*,*) 'roverw: R_inf/V_inf = ', R_inf/V_inf
 !  roverw = R_inf / V_inf
 END IF ! sobolev_approximation
ELSE
 write(*,*) 'roverw: velapprox = ', velapprox, ' is not a valid choice'
 STOP
END IF

END FUNCTION
