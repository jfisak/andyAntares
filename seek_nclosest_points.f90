! finds n_closest points 
!
! INPUT: cur_pos, DBLE(const_dimofspace), the current poisition 
!        cur_points, INT(cur_n_points): input indeces
!        n_closest, INT: number of the closest points we are looking for
!        cur_n_points, INT: size of the array cur_points
! OUTPUT: interp_dist: INT(ALLOC(:,:)): arrays containing index and distance
SUBROUTINE seek_nclosest_points(cur_pos, cur_points, n_closest, cur_n_points, closest_points)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_pos

DOUBLE PRECISION, ALLOCATABLE           :: pom(:,:)
DOUBLE PRECISION, DIMENSION(n_closest, 2) :: interp_dist
INTEGER, DIMENSION(n_closest)           :: closest_points

INTEGER, PARAMETER                      :: grid_A = 1, grid_B = 2
INTEGER                                 :: n_closest
INTEGER, PARAMETER                      :: ind_dist = 1, ind_index = 2
DOUBLE PRECISION, DIMENSION(3)          :: cur_mg_pos
INTEGER                                 :: cur_n_points

INTEGER                                 :: cur_nearest_point, cur_nop
DOUBLE PRECISION                        :: dist
LOGICAL                                 :: seeking, novyBod, nahrada
INTEGER                                 :: cip, cur_index_pos
INTEGER, DIMENSION(cur_n_points)        :: cur_points
DOUBLE PRECISION                        :: dist_A, dist_B, dist_C
! INTEGER                                 :: cur_index_i

DOUBLE PRECISION, PARAMETER             :: large_number = 1.D90

INTEGER                                 :: cur_iti_mgi
INTEGER, DIMENSION(3)                                   :: n_coor, count_xyz
INTEGER, PARAMETER                                      :: coor_x = 1, coor_y = 2, coor_z = 3
DOUBLE PRECISION                                        :: cur_dist

INTEGER                                                 :: cur_index_I, cur_index_J, cur_itj_mgi, cur_index_K
DOUBLE PRECISION, DIMENSION(3)                          :: cur_saved_mg_pos_A, cur_saved_mg_pos_B
DOUBLE PRECISION, DIMENSION(3)                          :: vec_AB, vec_AC

INTEGER, DIMENSION(3)                                   :: f_indexy, index_delete, mgi_indexy
INTEGER                                                 :: cur_del_index, cur_mgi

LOGICAL                                                 :: procout=.true.

! write(*,*) 'seek_nclosest_points: cur_points = ', cur_points
!_______________________________________________________________
!            SEEKING FOR THE n_closest CLOSEST POINTS
!_______________________________________________________________
 ! interp_dist(I,J)
 ! interp_dist(:,1) -- distance
 ! interp_dist(:,2) -- index of the point
interp_dist(:,ind_dist) = large_number
interp_dist(:,ind_index) = 0.D0
seeking = .TRUE.
! going through the array interp_dist and sorting from the smallest distances to largest distances
DO cur_nearest_point = 1, cur_n_points ! #L00
 IF(procout) write(*,*) 'seek_nclosest_points: _______________________________________________________________________'
 IF(procout) write(*,*) 'seek_nclosest_points: cur_point = ', cur_nearest_point, '/', cur_n_points, ' dist = ', dist
 cur_nop = cur_points(cur_nearest_point)
 IF(model_type == 2) THEN
  IF(inputmodel == 1) THEN
   cur_mg_pos(ind_x) = model_grid(cur_nop)%rwind
   cur_mg_pos(ind_y) = model_grid(cur_nop)%angle
   cur_mg_pos(ind_z) = 0.D0
  END IF
 ELSE IF(model_type == 3) THEN
  cur_mg_pos = model_grid(cur_nop)%vec_pos
 END IF
 ! write(45,*) cur_mg_pos
 ! distance of a current modGrid point and the 
 dist = sqrt((cur_pos(ind_x) - cur_mg_pos(ind_x))**2 + (cur_pos(ind_y) - cur_mg_pos(ind_y))**2 + &
                (cur_pos(ind_z) - cur_mg_pos(ind_z))**2)

 ! the current distance in the array is larger than the current distance, we should move all values
 ! one level below and save the new distance to the current index
 !write(*,*) 'seek_nclosest_points: intpdist = ', interp_dist(cur_index_pos, ind_dist), ' dist = ', dist
 seeking = .TRUE.
 ! testing of the coordinates
 n_coor(:) = 0
 novyBod = .true.
 nahrada = .false.
 index_delete(:) = 0
 f_indexy(:) = 0
 mgi_indexy(:) = 0
 ! going value by value and counting the number of points with very same coordinate
 ! the number larger than two is useless, therefore we skip this point
 write(*,*) 'seek_nclosest_points: cur_mg_pos = ', cur_mg_pos

 ! if the current point is not closer than the most far (last index) point in the interp_dist,
 ! we do not have to do anything and go to the following point
 IF(interp_dist(n_closest, ind_dist) < dist) THEN
  IF(procout) write(*,*) 'seek_nclosest_points: the point ', cur_nearest_point, ' is too far away'
  CYCLE
 END IF

 ! in this case, we have to go through the saved list o points to see, whether we can save the new point or not
 ! IF(n_coor(coor_x) >= 2 .or. n_coor(coor_y) >= 2 .or. n_coor(coor_z) >= 2) THEN
 count_xyz(:) = 0
 IF(procout) write(*,*) 'seek_nclosest_points: seeking for the points on the same line'
 cur_del_index = 1
 DO cur_index_I = 1, n_closest
  cur_iti_mgi = INT(interp_dist(cur_index_I, ind_index))
  DO cur_index_J = cur_index_I + 1, n_closest
   cur_itj_mgi = INT(interp_dist(cur_index_J, ind_index))
   IF(INT(cur_iti_mgi) /= 0 .and. INT(cur_itj_mgi) /= 0) THEN
    IF(procout) write(*,*) 'seek_nclosest_points: index = ', interp_dist(cur_index_I, ind_index), &
     interp_dist(cur_index_J, ind_index)
    IF(model_type == 2) THEN
     IF(inputmodel ==1) THEN
      cur_saved_mg_pos_A(ind_x) = model_grid(cur_iti_mgi)%rwind
      cur_saved_mg_pos_A(ind_y) = model_grid(cur_iti_mgi)%angle
      cur_saved_mg_pos_A(ind_z) = 0.D0
      cur_saved_mg_pos_B(ind_x) = model_grid(cur_itj_mgi)%rwind
      cur_saved_mg_pos_B(ind_y) = model_grid(cur_itj_mgi)%angle
      cur_saved_mg_pos_B(ind_z) = 0.D0
     END IF
    ELSE IF(model_type == 3) THEN
     cur_saved_mg_pos_A = model_grid(cur_iti_mgi)%vec_pos
     cur_saved_mg_pos_B = model_grid(cur_itj_mgi)%vec_pos
    END IF

    vec_AB = (cur_saved_mg_pos_A - cur_saved_mg_pos_B)/NORM2(cur_saved_mg_pos_A - cur_saved_mg_pos_B)
    vec_AC = (cur_saved_mg_pos_A - cur_mg_pos)/NORM2(cur_saved_mg_pos_A - cur_mg_pos)
    IF(procout) write(*,*) 'seek_nclosest_points: vektory vec_AB a vec_AC'
    IF(procout) write(*,*) 'vec_AB = ', vec_AB
    IF(procout) write(*,*) 'vec_AC = ', vec_AC
    ! testing if a new point is on the same line as the other points already saved into interp_dist
    IF(((vec_AB(1) == vec_AC(1)) .and. (vec_AB(2) == vec_AC(2)) .and. (vec_AB(3) == vec_AC(3))) .or. &
     & ((vec_AB(1) == -vec_AC(1)) .and. (vec_AB(2) == -vec_AC(2)) .and. (vec_AB(3) == -vec_AC(3)))) THEN 
     ! zatím to risknu a vezmu ten druhý index, který by měl být index vzdálenějšího bodu
     f_indexy(cur_del_index) = cur_index_J
     mgi_indexy(cur_del_index) = cur_itj_mgi
     cur_del_index = cur_del_index + 1
     IF(procout) write(*,*) 'seek_nclosest_points: for delete = ', cur_index_J
     novyBod = .FALSE.
    END IF ! vec1 = +- vec2
   END IF ! interp_dist(mod_index) /= 0
  END DO
 END DO ! cur_index_I = 1, n_closest
 
 IF(.not. novyBod) THEN
  ! looking for the less distant point
  IF(f_indexy(1) > 0) THEN
   dist_A = interp_dist(f_indexy(1), ind_dist)
  ELSE
   dist_A = large_number
  END IF
  IF(f_indexy(2) > 0) THEN
   dist_B = interp_dist(f_indexy(2), ind_dist)
  ELSE
   dist_B = large_number
  END IF
  IF(f_indexy(3) > 0) THEN
   dist_C = interp_dist(f_indexy(3), ind_dist)
  ELSE
   dist_C = large_number
  END IF

  ! all more distant points should be replaced with the current one
  IF(dist < dist_A .and. dist < dist_B .and. dist < dist_C) THEN
   IF(f_indexy(1) > 0) index_delete(1) = f_indexy(1)
   IF(f_indexy(2) > 0) index_delete(2) = f_indexy(2)
   IF(f_indexy(3) > 0) index_delete(3) = f_indexy(3)
   IF(procout) write(*,*) 'seek_nclosest_points: index_delete = ', index_delete(:)
   nahrada = .true.
  ELSE
   nahrada = .false.
  END IF
  IF(procout) write(*,*) 'seek_nclosest_points: index to delete: ', index_delete(:)

  ! we now delete the old point and add a better point into the array
  ! deleting a point on a position index_delete
  DO cur_index_K = 1,3
   IF(index_delete(cur_index_K) > 0) THEN
    if(procout) write(*,*) 'seek_nclosest_points: smazani bodu ', index_delete(cur_index_K), ' mod_grid = ',&
     & interp_dist(index_delete(cur_index_K), ind_index)
    DO cur_index_I = 1, n_closest
     cur_mgi = INT(interp_dist(cur_index_I, ind_index))
     IF(cur_mgi == mgi_indexy(cur_index_K)) THEN
      IF(procout) write(*,*) 'seek_nclosest_points: smazani indexu: ', cur_mgi
      DO cur_index_J = cur_index_I, n_closest
       IF(cur_index_J < n_closest) THEN
        interp_dist(cur_index_J, ind_dist) = interp_dist(cur_index_J + 1, ind_dist)
        interp_dist(cur_index_J, ind_index) = interp_dist(cur_index_J + 1, ind_index)
       ELSE IF(cur_index_J == n_closest) THEN
        interp_dist(cur_index_J, ind_dist) = large_number
        interp_dist(cur_index_J, ind_index) = 0
       END IF
      END DO
     END IF
    END DO
   END IF ! index_delete > 0
  END DO
 END IF ! not novyBod
  
 ! we can now put the new point into the list
 ! without erassing any other point
 write(*,*) 'seek_nclosest_points: novyBod = ', novyBod
 IF(novyBod .or. nahrada) THEN
  pom = interp_dist
  write(*,*) 'seek_nclosest_points: cur_index_pos = ', cur_index_pos, ' n_closest = ', n_closest
  DO cur_index_pos = 1, n_closest
   write(*,*) 'seek_nclosest_points: cip = ', cip, ' cur_index_pos = ', cur_index_pos
   cur_dist = interp_dist(cur_index_pos, ind_dist)
   write(*,*) 'seek_nclosest_points: cur_dist = ', cur_dist, ' dist = ', dist
   IF(dist < cur_dist) THEN
    DO cip = cur_index_pos, n_closest
     IF(cip == cur_index_pos) THEN
      interp_dist(cip, ind_dist) = dist
      interp_dist(cip, ind_index) = cur_nop
     ELSE
      interp_dist(cip, ind_dist) = pom(cip - 1, ind_dist)
      interp_dist(cip, ind_index) = pom(cip - 1, ind_index)
     END IF
     seeking = .FALSE.
    END DO
    EXIT
   END IF
  END DO
 END IF ! novyBod
 write(*,*) 'seek_nclosest_points: interp_dist() = ', INT(interp_dist(:,ind_index))
END DO ! #L00 loop over all possible points

closest_points = INT(interp_dist(:,ind_index))

write(*,*) 'seek_nclosest_points: interp_dist = ', interp_dist(:,ind_index)
write(*,*) 'seek_nclosest_points: interp_dist = ', interp_dist(:,ind_dist)

END SUBROUTINE seek_nclosest_points
