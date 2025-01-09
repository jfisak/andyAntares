! this sbr calculates a velocity vector via trilinear interpolation for the
! position of the packet
!
! INPUT: pack_index
! 
! OUTPUT: vel_vec
SUBROUTINE vel_discrete_points(pack_index, vel_vec)

USE types
USE constants
USE dummypacket
IMPLICIT NONE

INTEGER                                 :: pack_index, dummypack_index
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: vel_vec

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: act_corner, act_width, act_pos
INTEGER                                 :: act_cell
INTEGER                                 :: act_mgi
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: act_center, act_vel
DOUBLE PRECISION                        :: act_vel_norm

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: rel_pos

INTEGER, DIMENSION(8)                   :: velgridcells

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cur_pos1, cur_vel1, cur_pos2, cur_vel2
INTEGER                                 :: cur_cell1, cur_cell2, cur_mgi1, cur_mgi2

INTEGER                                 :: I, J
! DOUBLE PRECISION, DIMENSION(3,8)        :: c_point, c_pos
DOUBLE PRECISION, DIMENSION(const_dimofspace,4)        :: e_point, e_pos 
DOUBLE PRECISION, DIMENSION(const_dimofspace,2)        :: w_point, w_pos
! DOUBLE PRECISION, DIMENSION(3)          :: point_width

LOGICAL                                 :: incellmode

DOUBLE PRECISION, DIMENSION(8,const_dimofspace)        :: cube_pos

! INTEGER                                 :: cur_propGrid_cell


! write(*,*) 'vel_discrete_points: pack_index = ', pack_index, ' dim(package) = ', SIZE(package)
IF(pack_index <= SIZE(package)) THEN
 act_cell = package(pack_index)%cell_numb
 act_pos = package(pack_index)%pos
ELSE IF(pack_index > SIZE(package)) THEN
 dummypack_index = pack_index - SIZE(package)
 act_cell = dummypackage(dummypack_index)%cell_numb
 act_pos = dummypackage(dummypack_index)%pos
 ! CALL find_dyn_cell1(act_pos, cur_propGrid_cell)
 ! write(*,*) 'vel_discrete_points: act_cell = ', act_cell, ' cur_propGrid_cell = ', cur_propGrid_cell
END IF
act_mgi = dyn_cell(act_cell)%model_index
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width
act_center = act_corner + act_width/2.0

act_vel = model_grid(act_mgi)%vec_vel
act_vel_norm = model_grid(act_mgi)%vel

! a relative position in respect to the propCell center
rel_pos = act_pos - act_center


! cell neighbour numbers
IF(dyngrid == 0) THEN
 CALL oct_neighbors(pack_index, rel_pos, velgridcells, incellmode)
ELSE IF(dyngrid > 0) THEN
 CALL oct_virtcube(pack_index, rel_pos, cube_pos, incellmode)
END IF
! a special case when some neighbor cells do not exist, because we are close bound to the compuational domain
IF(incellmode) THEN ! incellmode
 CALL velo_vector(act_pos, act_mgi, vel_vec)
! trilinear interpolation
! z-direction -- four points
! IF(pack_index == 1) THEN
ELSE ! incellmode
 DO I = 1,4
  IF(dyngrid == 0) THEN
   cur_cell1 = velgridcells(2*I -1)
   cur_mgi1 = dyn_cell(cur_cell1)%model_index
   cur_pos1 = dyn_cell(cur_cell1)%corner + dyn_cell(cur_cell1)%width/2.0
   
   cur_cell2 = velgridcells(2*I)
   cur_mgi2 = dyn_cell(cur_cell2)%model_index
   cur_pos2 = dyn_cell(cur_cell2)%corner + dyn_cell(cur_cell2)%width/2.0
  ELSE IF (dyngrid > 0) THEN
   cur_pos1 = cube_pos(2*I - 1, :)
   CALL find_dyn_cell1(cur_pos1, cur_cell1)
   cur_mgi1 = dyn_cell(cur_cell1)%model_index

   cur_pos2 = cube_pos(2*I, :)
   CALL find_dyn_cell1(cur_pos2, cur_cell2)
   cur_mgi2 = dyn_cell(cur_cell2)%model_index
  END IF
  CALL velo_vector(cur_pos1, cur_mgi1, cur_vel1)

  CALL velo_vector(cur_pos2, cur_mgi2, cur_vel2)

  if(cur_vel1(ind_x) == -2.0 .and. cur_vel1(ind_y) == -3.0 .and. cur_vel1(ind_z) == -5.0) then
   cur_vel1 = act_vel
  end if
  if(cur_vel2(ind_x) == -2.0 .and. cur_vel2(ind_y) == -3.0 .and. cur_vel2(ind_z) == -5.0) then
   cur_vel2 = act_vel
  end if

  IF(cur_pos1(ind_z) > cur_pos2(ind_z) .and. (act_pos(ind_z) < cur_pos2(ind_z) .or. act_pos(ind_z) > cur_pos1(ind_z))) THEN
   DO J = 1,8
    write(29,*) cube_pos(I,:)
   END DO
  ELSE IF (cur_pos1(ind_z) < cur_pos2(ind_z) .and. (act_pos(ind_z) > cur_pos2(ind_z) .or. act_pos(ind_z) < cur_pos1(ind_z))) THEN
   DO J = 1,8
    write(29,*) cube_pos(I,:)
   END DO
  END IF
   

 ! IF(cur_pos1(3) == cur_pos2(3)) THEN
 !  DO J = 1,8
 !   write(29,*) cube_pos(I,:)
 !  END DO
 ! END IF
 
  CALL lin_interpolation(act_pos(ind_z), cur_vel1, cur_pos1(ind_z), cur_vel2, cur_pos2(ind_z), e_point(:,I))
  e_pos(:,I) = (/ cur_pos1(ind_x), cur_pos1(ind_y), act_pos(ind_z)/)
 END DO
 
 
 ! y-direction -- two points
 DO I = 1,2
  cur_vel1 = e_point(:, 2*I - 1)
  cur_pos1 = e_pos(:, 2*I - 1)
 
  cur_vel2 = e_point(:,2*I)
  cur_pos2 = e_pos(:,2*I)
 
  CALL lin_interpolation(act_pos(2), e_point(:,2*I-1), cur_pos1(2), e_point(:,2*I), cur_pos2(2), w_point(:,I))
  w_pos(:,I) = (/ cur_pos1(1), act_pos(2), e_pos(3,I) /)
 END DO
 
 ! x-direction -- one last point
  cur_vel1 = w_point(:,1)
  cur_pos1 = w_pos(:,1)
 
  cur_vel2 = w_point(:,2)
  cur_pos2 = w_pos(:,2)
  
  CALL lin_interpolation(act_pos(1), w_point(:,1), cur_pos1(1), w_point(:,2), cur_pos2(1), vel_vec)
END IF

! write(49,*) norm2(act_pos), norm2(vel_vec)



END SUBROUTINE vel_discrete_points
