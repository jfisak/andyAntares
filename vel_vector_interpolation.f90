! this sbr calculates a velocity vector for a selected point
! usually a centre of a propGrid cell based on the input velocity
! vectors
!
! input: positon -- the position of the interpolation
!        mgi_indexes -- the indexes of the closest model grids
!        n_clo_mgi -- the number of the closest mgi points
! output: velocity -- velocity vector
SUBROUTINE vel_vector_interpolation(positon, mgi_indexes, n_clo_mgi, velocity)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: positon, velocity
INTEGER                                                 :: n_clo_mgi, n_clo_mgi_half
DOUBLE PRECISION, DIMENSION(n_clo_mgi)                  :: mgi_indexes

DOUBLE PRECISION, DIMENSION(3)                          :: summ_r, cur_pos, centre
INTEGER                                                 :: cur_I

! sorting
INTEGER                                                 :: cur_index, cur_iter, cur_mgi, cur_mgi_iter
DOUBLE PRECISION                                        :: dummy_var, dummy

INTEGER, PARAMETER                                      :: coor_x = 1, coor_y = 2, coor_z = 3
INTEGER                                                 :: cur_coordinate, cur_mgi1, cur_mgi2, cur_pair
DOUBLE PRECISION, DIMENSION(3)                          :: cur_pos1, cur_pos2, cur_vel1, cur_vel2, point_pos
DOUBLE PRECISION, ALLOCATABLE                           :: vel_8(:,:), pos_8(:,:), vel_4(:,:), pos_4(:,:)

INTEGER                                                 :: cur_mgiind

! write(*,*) 'vel_vector_interpolation: mgi_indexes = ', mgi_indexes(:)

! the ``centre of the mass'' of the points
summ_r(:) = 0.D0
DO cur_I = 1, n_clo_mgi
 cur_mgiind = INT(mgi_indexes(cur_I))
 cur_pos = model_grid(cur_mgiind)%vec_pos
 summ_r(1) = summ_r(1) + cur_pos(1)
 summ_r(2) = summ_r(2) + cur_pos(2)
 summ_r(3) = summ_r(3) + cur_pos(3)
END DO

centre = summ_r / n_clo_mgi

! write(45,*) centre/R_star

! sorting points according the z-coordinate

DO cur_index = 2, n_clo_mgi
 cur_iter = cur_index - 1

 cur_mgi = INT(mgi_indexes(cur_index))
 cur_mgi_iter = INT(mgi_indexes(cur_iter))

 dummy_var = model_grid(cur_mgi)%vec_pos(3)

 DO WHILE(cur_iter >= 1)
  IF(model_grid(cur_mgi_iter)%vec_pos(3) > dummy_var) THEN
   dummy = INT(mgi_indexes(cur_iter + 1))
   mgi_indexes(cur_iter + 1) = mgi_indexes(cur_iter)
   mgi_indexes(cur_iter) = dummy
  END IF
  cur_iter = cur_iter - 1
 END DO
END DO


! now we have n/2 + n/2 indeces which we sort according the y coordinate
n_clo_mgi_half = n_clo_mgi / 2

! the first half
DO cur_index = 2, n_clo_mgi/2
 cur_iter = cur_index - 1

 cur_mgi = INT(mgi_indexes(cur_index))
 cur_mgi_iter = INT(mgi_indexes(cur_iter))

 dummy_var = model_grid(cur_mgi)%vec_pos(2)

 DO WHILE(cur_iter >= 1)
  IF(model_grid(cur_mgi_iter)%vec_pos(2) > dummy_var) THEN
   dummy = mgi_indexes(cur_iter + 1)
   mgi_indexes(cur_iter + 1) = mgi_indexes(cur_iter)
   mgi_indexes(cur_iter) = dummy
  END IF
  cur_iter = cur_iter - 1
 END DO
END DO

! the second half
DO cur_index = n_clo_mgi_half + 2, n_clo_mgi
 cur_iter = cur_index - 1

 cur_mgi = INT(mgi_indexes(cur_index))
 cur_mgi_iter = INT(mgi_indexes(cur_iter))

 dummy_var = model_grid(cur_mgi)%vec_pos(2)

 DO WHILE(cur_iter >= n_clo_mgi_half + 1)
  IF(model_grid(cur_mgi_iter)%vec_pos(2) > dummy_var) THEN
   dummy = mgi_indexes(cur_iter + 1)
   mgi_indexes(cur_iter + 1) = mgi_indexes(cur_iter)
   mgi_indexes(cur_iter) = dummy
  END IF
  cur_iter = cur_iter - 1
 END DO
END DO

DO cur_index = 1, n_clo_mgi
 cur_mgi = INT(mgi_indexes(cur_index))
 ! write(*,*) 'vel_vector_interpolation: pos = ', model_grid(cur_mgi)%vec_pos
END DO

ALLOCATE(vel_8(4,3), pos_8(4,3), vel_4(2,3), pos_4(2,3))

! trilinear interpolation
! eight points
cur_coordinate = coor_x
DO cur_pair = 1, 4
 cur_mgi1 = INT(mgi_indexes(2 * cur_pair - 1))
 cur_vel1 = model_grid(cur_mgi1)%vec_vel
 cur_pos1 = model_grid(cur_mgi1)%vec_pos

 cur_mgi2 = INT(mgi_indexes(2 * cur_pair))
 cur_vel2 = model_grid(cur_mgi2)%vec_vel
 cur_pos2 = model_grid(cur_mgi2)%vec_pos

 CALL lin_interpolation_3D(positon, cur_pos1, cur_pos2, cur_vel1, cur_vel2, cur_coordinate, velocity)
 ! write(43,*) cur_pos1, cur_vel1
 ! write(43,*) cur_pos2, cur_vel2
 point_pos = cur_pos1 + (cur_pos2 - cur_pos1)* &
   & (positon(cur_coordinate) - cur_pos1(cur_coordinate))/(cur_pos2(cur_coordinate) - cur_pos1(cur_coordinate))
 IF(isnan(point_pos(1)) .or. isnan(point_pos(2)) .or. isnan(point_pos(3))) THEN
  write(*,*) 'vel_vector_interpolation: pos2(x) = ', cur_pos2(cur_coordinate), ' pos1(x) = ', cur_pos1(cur_coordinate)
  DO  cur_I = 1, n_clo_mgi
   cur_mgi = INT(mgi_indexes(cur_I))
   cur_vel1 = model_grid(cur_mgi)%vec_vel
   cur_pos1 = model_grid(cur_mgi)%vec_pos
   write(43,*) cur_pos1, cur_vel1
  END DO
  write(44,*) positon
  STOP 'vel_vector_interpolation: NaNs'
 END IF
 ! saving calculated quantities into an array
 vel_8(cur_pair,:) = velocity(:)
 pos_8(cur_pair,:) = point_pos(:)
END DO

! write(45,*) positon
! trilinear interpolation
! four points
cur_coordinate = coor_y
DO cur_pair = 1,2
 cur_vel1 = vel_8(2*cur_pair - 1,:)
 cur_pos1 = pos_8(2*cur_pair - 1,:)

 cur_vel2 = vel_8(2*cur_pair,:)
 cur_pos2 = pos_8(2*cur_pair,:)

 write(*,*) 'vel_vector_interpolation: cur_pos1 = ', cur_pos1, ' cur_pos2 = ', cur_pos2
 write(*,*) 'vel_vector_interpolation: cur_vel1 = ', cur_vel1, ' cur_vel2 = ', cur_vel2

 CALL lin_interpolation_3D(positon, cur_pos1, cur_pos2, cur_vel1, cur_vel2, cur_coordinate, velocity)
 point_pos = cur_pos1 + (cur_pos2 - cur_pos1)* &
   & (positon(cur_coordinate) - cur_pos1(cur_coordinate))/(cur_pos2(cur_coordinate) - cur_pos1(cur_coordinate))
 vel_4(cur_pair,:) = velocity(:)
 pos_4(cur_pair,:) = point_pos(:)
END DO

! the last interpolation
 cur_pair = 1
 cur_coordinate = coor_z
 cur_vel1 = vel_4(2*cur_pair - 1,:)
 cur_pos1 = pos_4(2*cur_pair - 1,:)

 cur_vel2 = vel_4(2*cur_pair,:)
 cur_pos2 = pos_4(2*cur_pair,:)

 CALL lin_interpolation_3D(positon, cur_pos1, cur_pos2, cur_vel1, cur_vel2, cur_coordinate, velocity)



END SUBROUTINE vel_vector_interpolation
