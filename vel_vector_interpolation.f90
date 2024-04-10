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


! the ``centre of the mass'' of the points
summ_r(:) = 0.D0
DO cur_I = 1, n_clo_mgi
 cur_pos = model_grid(mgi_indexes(cur_I))%vec_pos
 summ_r(1) = summ_r(1) + cur_pos(1)
 summ_r(2) = summ_r(2) + cur_pos(2)
 summ_r(3) = summ_r(3) + cur_pos(3)
END DO

centre = summ_r / n_clo_mgi

write(45,*) centre/R_star

! sorting points according the z-coordinate

DO cur_index = 2, n_clo_mgi
 cur_iter = cur_index - 1

 cur_mgi = mgi_indexes(cur_index)
 cur_mgi_iter = mgi_indexes(cur_iter)

 dummy_var = model_grid(cur_mgi)%vec_pos(3)

 DO WHILE(cur_iter >= 1)
  IF(model_grid(cur_mgi_iter)%vec_pos(3) > dummy_var) THEN
   dummy = mgi_indexes(cur_iter + 1)
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

 cur_mgi = mgi_indexes(cur_index)
 cur_mgi_iter = mgi_indexes(cur_iter)

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

 cur_mgi = mgi_indexes(cur_index)
 cur_mgi_iter = mgi_indexes(cur_iter)

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
 cur_mgi = mgi_indexes(cur_index)
 write(*,*) 'vel_vector_interpolation: pos = ', model_grid(cur_mgi)%vec_pos
END DO


! trilinear interpolation
! eight points
cur_coordinate = coor_x
DO cur_pair = 1, 4
 cur_mgi1 = mgi_indexes(2 * cur_pair - 1)
 cur_vel1 = model_grid(cur_mgi1)%vec_vel
 cur_pos1 = model_grid(cur_mgi1)%vec_pos

 cur_mgi2 = mgi_indexes(2 * cur_pair)
 cur_vel2 = model_grid(cur_mgi2)%vec_vel
 cur_pos2 = model_grid(cur_mgi2)%vec_pos

 CALL lin_interpolation_3D()

END DO







END SUBROUTINE vel_vector_interpolation
