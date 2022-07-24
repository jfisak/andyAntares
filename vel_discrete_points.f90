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
DOUBLE PRECISION, DIMENSION(3,8)        :: c_point, c_pos
DOUBLE PRECISION, DIMENSION(3,4)        :: e_point, e_pos 
DOUBLE PRECISION, DIMENSION(3,2)        :: w_point, w_pos
DOUBLE PRECISION, DIMENSION(3)          :: point_width

INTEGER                                 :: pomocna_bunka
LOGICAL                                 :: incellmode

INTEGER                                 :: dummypackage

dummypackage = SIZE(package)


testPacket = SIZE(package) - 1
package(testPacket) = package(pack_index)
act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
act_mgi = get_package_model_index(pack_index)
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width
act_center = act_corner + act_width/2.0

act_vel_norm = model_grid(act_mgi)%vel
! write(*,*) 'vel_discrete_points: act_vel_norm = ', act_vel_norm/V_inf

! a relative position in respect to the propCell center
rel_pos = act_pos - act_center
! write(*,*) 'vel_discrete_points: rel_pos = ', rel_pos(:)

! DO I = 1,3
!  IF((act_pos(I) < act_corner(I) .OR. act_pos(I) > act_corner(I) + act_width(I) ) .and. pack_index /= dummypackage) THEN
!   CALL find_dyn_cell1(act_pos, pomocna_bunka)
!   write(*,*) 'writing into the file fort.4'
!   write(*,*) 'vel_discrete_points: act_cell = ', act_cell, ' skutecna bunka = ', pomocna_bunka
!   ! write(4,*) act_pos, act_corner, act_width, get_package_model_index(pack_index)
!   ! write(4,*) act_pos, dyn_cell(pomocna_bunka)%corner, dyn_cell(pomocna_bunka)%width, get_package_model_index(pack_index)
!   write(*,*) 'vel_discrete_points: pack_index = ', pack_index
!   STOP 'vel_discrete_points: packet is not located inside the propagation cell'
!  END IF
! END DO



! cell neighbour numbers
CALL oct_neighbors(pack_index, rel_pos, velgridcells, incellmode)
!  DO I = 1,8
!   write(25,*) dyn_cell(velgridcells(I))%corner + dyn_cell(velgridcells(I))%width/2.0
!  END DO
!  write(25,*) act_pos
! write(*,*) 'vel_discrete_points: incellmode = ', incellmode
! a special case when some neighbor cells do not exist, because we are close bound to the compuational domain
! we will choose initial points on the bound of the current propagation grid instead
IF(incellmode) THEN ! incellmode
 point_width = act_width
 IF(rel_pos(3) < 0) THEN
  point_width(3) = -act_width(3)
 END IF
 IF(rel_pos(2) < 0) THEN
  point_width(2) = -act_width(2)
 END IF
 IF(rel_pos(1) < 0) THEN
  point_width(1) = -act_width(1)
 END IF
 ! points for calculation
 DO I = 1,8
  c_pos(:,I) = act_center
 END DO
 ! +z
 c_pos(3,2) = act_center(3) + point_width(3)
 ! +y
 c_pos(2,3) = act_center(2) + point_width(2)
 ! +y+z
 c_pos(2,4) = act_center(2) + point_width(2)
 c_pos(3,4) = act_center(3) + point_width(3)
 ! +x
 c_pos(1,5) = act_center(1) + point_width(1)
 ! +x+z
 c_pos(1,6) = act_center(1) + point_width(1)
 c_pos(3,6) = act_center(3) + point_width(3)
 ! +x+y
 c_pos(1,7) = act_center(1) + point_width(1)
 c_pos(2,7) = act_center(2) + point_width(2)
 ! +x+y+z
 c_pos(1,8) = act_center(1) + point_width(1)
 c_pos(2,8) = act_center(2) + point_width(2)
 c_pos(3,8) = act_center(3) + point_width(3)

  ! write(*,*) 'vel_discrete_points: ', 'point_width = ', point_width, ' c_pos = ', c_pos
  ! DO I = 1,8
  !  write(28,*) c_pos(:,I)
  ! END DO
 ! write(27,*) act_pos
 ! zero approximation
 ! now we expect all vectors to be equal to zero, except the center vector
 IF(model_type == 1) c_point(:,1) = act_vel_norm * c_pos(:,1)/norm2(c_pos(:,1))
 DO I = 2,8
  c_point(:,I) = (/0.0, 0.0, 0.0 /)
 END DO

 
 DO I = 1,4
  cur_pos1 = c_pos(:,2*I-1)
  cur_vel1 = c_point(:,2*I-1)
  
  cur_pos2 = c_pos(:,2*I)
  cur_vel2 = c_point(:,2*I)
 
  ! write(*,*) 'vel_discrete_points: cur_vel1 = ', cur_vel1, ' cur_vel2 = ', cur_vel2
  ! write(*,*) 'vel_discrete_points: e pack_index = ', pack_index, ' cur_pos1 = ', cur_pos1(3), ' cur_pos2 = ', cur_pos2(3)
  ! CALL lin_interpolation(act_pos(3), cur_vel1, cur_pos1(3), cur_vel2, cur_pos2(3), e_point(:,I), act_pos, &
  CALL lin_interpolation(act_pos(3), cur_vel1, cur_pos1(3), cur_vel2, cur_pos2(3))
  e_pos(:,I) = (/ cur_pos1(1), cur_pos1(2), act_pos(3)/)
  ! write(*,*) 'vel_discrete_points: e_point = ', e_point(:,I)
  ! write(26,*) e_pos(:,I), e_point(:,I)
  ! write(*,*) 'vel_discrete_points: e/c = ', norm2(e_point(:,I))/light_speed
 END DO



 DO I = 1,2
  cur_vel1 = e_point(:,2*I-1)
  cur_pos1 = e_pos(:,2*I - 1)
 
  cur_vel2 = e_point(:,2*I)
  cur_pos2 = e_pos(:,2*I)
 
  ! write(*,*) 'vel_discrete_points: w cur_pos1 = ', cur_pos1(2), ' cur_pos2 = ', cur_pos2(2)
  CALL lin_interpolation(act_pos(2), e_point(:,2*I-1), cur_pos1(2), e_point(:,2*I), cur_pos2(2))
  w_pos(:,I) = (/ cur_pos1(1), act_pos(2), e_pos(3,I) /)
  ! write(25,*) w_pos(:,I), w_point(:,I)
  ! write(*,*) 'vel_discrete_points: w/c = ', norm2(w_point(:,I))/light_speed
 END DO
 
 ! x-direction -- one last point
  cur_vel1 = w_point(:,1)
  cur_pos1 = w_pos(:,1)
 
  cur_vel2 = w_point(:,2)
  cur_pos2 = w_pos(:,2)
  
  ! write(*,*) 'vel_discrete_points: f cur_pos1 = ', cur_pos1(1), ' cur_pos2 = ', cur_pos2(1)
  CALL lin_interpolation(act_pos(1), w_point(:,1), cur_pos1(1), w_point(:,2), cur_pos2(1))
  ! write(*,*) 'vel_discrete_points: vel_vec = ', vel_vec
  ! write(*,*) 'vel_discrete_points: v = ', norm2(vel_vec)/light_speed
 ! STOP 'vel_discrete_points: testing'



! trilinear interpolation
! z-direction -- four points
! IF(pack_index == 1) THEN
ELSE ! incellmode
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
  ! write(*,*) 'vel_discrete_points: e cur_pos1 = ', cur_pos1, ' cur_pos2 = ', cur_pos2
  CALL lin_interpolation(act_pos(3), cur_vel1, cur_pos1(3), cur_vel2, cur_pos2(3))
  e_pos(:,I) = (/ cur_pos1(1), cur_pos1(2), act_pos(3 )/)
  ! write(*,*) 'vel_discrete_points: e_point = ', e_point(:,I)
  IF(pack_index == 1) write(28,*) e_pos(:,I), e_point(:,I)
  ! write(*,*) 'vel_discrete_points: e/c = ', norm2(e_point(:,I))/light_speed
 END DO
 
 
 ! y-direction -- two points
 DO I = 1,2
  cur_vel1 = e_point(:,2*I-1)
  cur_pos1 = e_pos(:,2*I - 1)
 
  cur_vel2 = e_point(:,2*I)
  cur_pos2 = e_pos(:,2*I)
 
  ! write(*,*) 'vel_discrete_points: w cur_pos1 = ', cur_pos1, ' cur_pos2 = ', cur_pos2

  CALL lin_interpolation(act_pos(2), e_point(:,2*I-1), cur_pos1(2), e_point(:,2*I), cur_pos2(2))
  w_pos(:,I) = (/ cur_pos1(1), act_pos(2), e_pos(3,I) /)
  ! write(27,*) w_pos(:,I), w_point(:,I)
  ! write(*,*) 'vel_discrete_points: w/c = ', norm2(w_point(:,I))/light_speed
 END DO
 
 ! x-direction -- one last point
  cur_vel1 = w_point(:,1)
  cur_pos1 = w_pos(:,1)
 
  cur_vel2 = w_point(:,2)
  cur_pos2 = w_pos(:,2)
  
  CALL lin_interpolation(act_pos(1), w_point(:,1), cur_pos1(1), w_point(:,2), cur_pos2(1))
END IF



END SUBROUTINE vel_discrete_points
