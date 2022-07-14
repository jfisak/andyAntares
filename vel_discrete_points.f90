SUBROUTINE vel_discrete_points(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

INTEGER                                 :: testpacket
DOUBLE PRECISION, DIMENSION(3)          :: act_corner, act_width, act_pos
INTEGER                                 :: act_cell
INTEGER                                 :: get_package_model_index, act_mgi
DOUBLE PRECISION, DIMENSION(3)          :: act_center
DOUBLE PRECISION                        :: act_vel_norm

DOUBLE PRECISION, DIMENSION(3)          :: rel_pos

INTEGER, DIMENSION(8)                   :: velgridcells

DOUBLE PRECISION, DIMENSION(3)          :: cur_pos1, cur_vel1, cur_pos2, cur_vel2
INTEGER                                 :: cur_cell1, cur_cell2, cur_mgi1, cur_mgi2

INTEGER                                 :: I
DOUBLE PRECISION, DIMENSION(3,4)        :: e_point, e_pos 
DOUBLE PRECISION, DIMENSION(3,2)        :: w_point, w_pos

testPacket = SIZE(package) - 1
package(testPacket) = package(pack_index)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
act_mgi = get_package_model_index(pack_index)
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width
act_center = act_corner + act_width/2.0

act_vel_norm = model_grid(act_mgi)%vel

! a relative position in respect to the propCell center
rel_pos = act_pos - act_center

! cell neighbour numbers
CALL oct_neighbors(pack_index, rel_pos, velgridcells)

! trilinear interpolation
! z-direction -- four points
DO I = 1,4
 cur_cell1 = velgridcells(2*I -1)
 cur_mgi1 = dyn_cell(cur_cell1)%model_index
 cur_pos1 = dyn_cell(cur_cell1)%corner + dyn_cell(cur_cell1)%width/2.0
 ! only for spherically symmetric case
 IF(model_type == 1) cur_vel1 = model_grid(cur_mgi1)%vel*cur_pos1/norm2(cur_pos1)
 
 cur_cell2 = velgridcells(2*I)
 cur_mgi2 = dyn_cell(cur_cell2)%model_index
 cur_pos2 = dyn_cell(cur_cell2)%corner + dyn_cell(cur_cell2)%width/2.0
 ! only for spherically symmetric case
 IF(model_type == 1) cur_vel2 = model_grid(cur_mgi2)%vel*cur_pos2/norm2(cur_pos2)

 ! write(*,*) 'vel_discrete_points: cur_vel1 = ', cur_vel1, ' cur_vel2 = ', cur_vel2
 CALL lin_interpolation(act_pos(3), cur_vel1, cur_pos1(3), cur_vel2, cur_pos2(3), e_point(:,I))
 e_pos(:,I) = (/ cur_pos1(1), cur_pos1(2), act_pos(3 )/)
 write(*,*) 'vel_discrete_points: e_point = ', e_point(:,I)
END DO

! y-direction -- two points
DO I = 1,2
 cur_vel1 = e_point(:,2*I-1)
 cur_pos1 = e_pos(:,2*I - 1)

 cur_vel2 = e_point(:,2*I)
 cur_pos2 = e_pos(:,2*I)

 CALL lin_interpolation(act_pos(2), cur_vel1, cur_pos1(2), cur_vel2, cur_pos2(2), w_point(:,I))
 w_pos(:,I) = (/ cur_pos1(1), act_pos(2), e_point(3,I) /)
END DO

! x-direction -- one last point
 cur_vel1 = w_point(:,1)
 cur_pos1 = w_pos(:,1)

 cur_vel2 = w_point(:,2)
 cur_pos2 = w_pos(:,2)
 
 CALL lin_interpolation(act_pos(1), cur_vel1, cur_pos1(1), cur_vel2, cur_pos2(1), vel_vec)

 write(*,*) 'vel_discrete_points: vel_vec = ', vel_vec





 STOP 'vel_discrete_points: testing'





END SUBROUTINE vel_discrete_points
