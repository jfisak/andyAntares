! this sbr interpolates the vector field of velocity from the modGrid into the propGrid
! 
! input: none
!
! output: saves new velocity vector for each* propGrid cell
!
!
SUBROUTINE vel_interpolation()

USE types
IMPLICIT NONE

INTEGER                                 :: N_vgrid_x, N_vgrid_y, N_vgrid_z
INTEGER                                 :: N_vgrid_cells_A, N_vgrid_cells_B
DOUBLE PRECISION                        :: w_vgrid_x, w_vgrid_y, w_vgrid_z
INTEGER                                 :: cur_point
INTEGER                                 :: n_A, n_B

INTEGER, PARAMETER                      :: min_incell = 10
INTEGER                                 :: n_in_cell

INTEGER, ALLOCATABLE                    :: n_points_A(:), n_points_B(:), &
                                         & indices_A(:), indices_B(:)
INTEGER, DIMENSION(n_modelgrid, 2)      :: vg_indexy_A, vg_indexy_B

DOUBLE PRECISION                        :: cur_x, cur_y, cur_z
INTEGER                                 :: cur_n_x_A, cur_n_y_A, cur_n_z_A
INTEGER                                 :: cur_n_x_B, cur_n_y_B, cur_n_z_B

INTEGER                                 :: cur_iter

INTEGER, DIMENSION(2)                   :: dummy_var_A, dummy_var_B, dummy_A, dummy_B
INTEGER                                 :: cur_ind_A, cur_ind_B, cur_vpg_cell

INTEGER                                 :: n_zeros

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cur_pg_pos, cur_mg_pos
INTEGER                                 :: cur_prop_cell

DOUBLE PRECISION, PARAMETER             :: large_number = 1.D90

INTEGER, ALLOCATABLE                    :: cur_points(:)
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cur_center_A, cur_center_B
INTEGER                                 :: cur_start_index, cur_end_index, cur_n_points
DOUBLE PRECISION                        :: dist_A, dist_B, dist_C

INTEGER, PARAMETER                      :: n_closest = 8
DOUBLE PRECISION, ALLOCATABLE           :: interp_dist(:,:), pom(:,:)

INTEGER                                 :: chosen_gridAB
INTEGER, PARAMETER                      :: grid_A = 1, grid_B = 2
INTEGER, PARAMETER                      :: ind_dist = 1, ind_index = 2

INTEGER                                 :: cur_nearest_point, cur_nop, cur_vg_index, cur_vg_point
DOUBLE PRECISION                        :: dist
LOGICAL                                 :: seeking, novyBod, nahrada
INTEGER                                                 :: cur_index, cip, cur_index_pos
! INTEGER                                 :: cur_index_i

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: vel_vector

INTEGER                                 :: cur_iti_mgi
INTEGER, DIMENSION(3)                                   :: n_coor, count_xyz
INTEGER, PARAMETER                                      :: coor_x = 1, coor_y = 2, coor_z = 3
DOUBLE PRECISION                                        :: cur_dist

INTEGER                                                 :: cur_index_I, cur_index_J, cur_itj_mgi, cur_index_K
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: cur_saved_mg_pos_A, cur_saved_mg_pos_B
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: vec_AB, vec_AC

INTEGER, DIMENSION(3)                                   :: f_indexy, index_delete, mgi_indexy
INTEGER                                                 :: cur_del_index, cur_mgi

LOGICAL                                                 :: procout=.true.
!_______________________________________________________________
!   #00             SET UP OF VIRTUAL GRIDS
!_______________________________________________________________
n_in_cell = FLOOR((n_modelgrid/min_incell)**(1.0/3.0))

! the size of grid should be at least equal to one
IF(n_in_cell < 1) n_in_cell = 1

N_vgrid_x = n_in_cell
N_vgrid_y = n_in_cell
N_vgrid_z = n_in_cell
N_vgrid_cells_A = N_vgrid_x * N_vgrid_y * N_vgrid_z
N_vgrid_cells_B = (N_vgrid_x - 1) * (N_vgrid_y - 1) * (N_vgrid_z - 1)
ALLOCATE(n_points_A(N_vgrid_cells_A), n_points_B(N_vgrid_cells_B))
ALLOCATE(indices_A(N_vgrid_cells_A), indices_B(N_vgrid_cells_B))
n_points_A(:) = 0
n_points_B(:) = 0
n_zeros = 0

w_vgrid_x = abs(xmax - xmin)/N_vgrid_x
w_vgrid_y = abs(ymax - ymin)/N_vgrid_y
w_vgrid_z = abs(zmax - zmin)/N_vgrid_z

!_______________________________________________________________
!   #01        CALCULATE VG INDEX OF MG CELLS
!_______________________________________________________________
! calculation of the virGrid index
DO cur_point = 1, n_modelgrid
 cur_x = model_grid(cur_point)%vec_pos(ind_x)
 cur_y = model_grid(cur_point)%vec_pos(ind_y)
 cur_z = model_grid(cur_point)%vec_pos(ind_z)
 cur_n_x_A = floor((cur_x - xmin)/w_vgrid_x) + 1
 cur_n_y_A = floor((cur_y - ymin)/w_vgrid_y) + 1
 cur_n_z_A = floor((cur_z - zmin)/w_vgrid_z) + 1
 cur_center_A = (/ w_vgrid_x * (cur_n_x_A + 5.D-1) + xmin, w_vgrid_y * (cur_n_y_A + 5.D-1) + ymin, &
                 & w_vgrid_z * (cur_n_z_A + 5.D-1) + zmin /)
 ! write(*,*) 'vel_interpolation: cur_center_A = ', cur_center_A

 ! n_A -- numerical index of VG cell
 n_A = cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1) + N_vgrid_x * N_vgrid_y * (cur_n_z_A - 1)
 ! write(*,*) 'vel_interpolation: n_x = ', cur_n_x_A, ' n_y = ', cur_n_y_A, ' n_z = ', cur_n_z_A
 ! write(*,*) 'vel_interpolation: n_A = ', n_A
 ! n_points -- number of points for the given cell
 n_points_A(n_A) = n_points_A(n_A) + 1
 ! vg_indexy -- list of indeces model grid --> VG index point
 vg_indexy_A(cur_point, 1) = n_A
 vg_indexy_A(cur_point, 2) = cur_point
 !!!!!!!!
 ! repete for the B grid
 IF(cur_x > xmin + w_vgrid_x/2.0 .and. cur_x < xmax - w_vgrid_x/2.0 .and.&
  & cur_y > ymin + w_vgrid_y/2.0 .and. cur_y < ymax - w_vgrid_y/2.0 .and. &
  & cur_z > zmin + w_vgrid_z/2.0 .and. cur_z < zmax - w_vgrid_z/2.0 ) THEN
  cur_n_x_B = floor((cur_x - xmin)/w_vgrid_x - 1.0/2.0) + 1
  cur_n_y_B = floor((cur_y - ymin)/w_vgrid_y - 1.0/2.0) + 1
  cur_n_z_B = floor((cur_z - zmin)/w_vgrid_z - 1.0/2.0) + 1
  n_B = cur_n_x_B + (N_vgrid_x - 1) * (cur_n_y_B - 1) + (N_vgrid_x - 1) * (N_vgrid_y - 1) * (cur_n_z_B - 1)
  n_points_B(n_B) = n_points_B(n_B) + 1
  vg_indexy_B(cur_point, 1) = n_B
  vg_indexy_B(cur_point, 2) = cur_point
 ELSE
  vg_indexy_B(cur_point, 1) = 0
  vg_indexy_B(cur_point, 2) = cur_point
  n_zeros = n_zeros + 1
 END IF
END DO

!_______________________________________________________________
!    #02            SORTING
!_______________________________________________________________
! sort the vg_indexy according to the VG index
DO cur_point = 2, n_modelgrid
 cur_iter = cur_point - 1

 dummy_var_A = vg_indexy_A(cur_point,:)

 DO WHILE(cur_iter >= 1)
  IF(vg_indexy_A(cur_iter,1) > dummy_var_A(1)) THEN
   ! A grid
   dummy_A = vg_indexy_A(cur_iter + 1,:)
   vg_indexy_A(cur_iter + 1,:) = vg_indexy_A(cur_iter,:)
   vg_indexy_A(cur_iter,:) = dummy_A
  END IF
  cur_iter = cur_iter - 1
 END DO
END DO

!_______________________________________________________________
! sort the vg_indexy according to the VG index
DO cur_point = 2, n_modelgrid
 cur_iter = cur_point - 1
 dummy_var_B = vg_indexy_B(cur_point,:)
 !_______________________________________________________________
 DO WHILE(cur_iter >= 1)
  IF(vg_indexy_B(cur_iter,1) > dummy_var_B(1)) THEN
   ! B grid
   dummy_B = vg_indexy_B(cur_iter + 1,:)
   vg_indexy_B(cur_iter + 1,:) = vg_indexy_B(cur_iter,:)
   vg_indexy_B(cur_iter,:) = dummy_B
  END IF
  cur_iter = cur_iter - 1
 END DO
END DO
!_______________________________________________________________

!_______________________________________________________________
! index array
!_______________________________________________________________

cur_ind_A = 0
cur_ind_B = n_zeros
!_______________________________________________________________
! create arrays with indeces pointing to an ordered list of modCell grids indeces
DO cur_vpg_cell = 1, N_vgrid_cells_A
 indices_A(cur_vpg_cell) = cur_ind_A + 1
 cur_ind_A = cur_ind_A + n_points_A(cur_vpg_cell)
END DO
!_______________________________________________________________
DO cur_vpg_cell = 1, N_vgrid_cells_B
 indices_B(cur_vpg_cell) = cur_ind_B + 1
 cur_ind_B = cur_ind_B + n_points_B(cur_vpg_cell)
END DO

ALLOCATE(interp_dist(n_closest,ind_index), pom(n_closest, ind_index))
! write(*,*) 'vel_interpolation: indices_B = ', indices_B
!_______________________________________________________________
!  #03              INTERPOLATING
!_______________________________________________________________

! going through the propGrid and finding a velocity vector using the trilinear interpolation
DO cur_prop_cell = 1, n_propgcells
 ! write(*,*) 'vel_interpolation: cur_prop_cell = ', cur_prop_cell
 !_______________________________________________________________
 !  #              CHOISE OF GRID A OR B
 !_______________________________________________________________
 ! the index in the AB grid
 ! cur_prop_cell = 1500
 cur_pg_pos = dyn_cell(cur_prop_cell)%corner + dyn_cell(cur_prop_cell)%width/2.D0
 cur_n_x_A = floor((cur_pg_pos(ind_x) - xmin)/w_vgrid_x) + 1
 cur_n_y_A = floor((cur_pg_pos(ind_y) - ymin)/w_vgrid_y) + 1
 cur_n_z_A = floor((cur_pg_pos(ind_z) - zmin)/w_vgrid_z) + 1
 n_A = cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1) + N_vgrid_x * N_vgrid_y * (cur_n_z_A - 1)
 cur_center_A = (/ w_vgrid_x * (cur_n_x_A + 5.D-1) + xmin, w_vgrid_y * (cur_n_y_A + 5.D-1) + ymin, &
                 & w_vgrid_z * (cur_n_z_A + 5.D-1) + zmin /)
 ! write(*,*) 'vel_interpolation: cur_center_A = ', cur_center_A
 ! write(*,*) 'vel_interpolation: N_vgrid_x = ', N_vgrid_x, ' N_vgrid_y = ', N_vgrid_y, ' N_vgrid_z = ', N_vgrid_z
 ! write(*,*) 'vel_interpolation: w_vgrid_x = ', w_vgrid_x/R_star, ' w_vgrid_y = ', w_vgrid_y/R_star, &
 !               & ' w_vgrid_z = ', w_vgrid_z/R_star
 ! write(*,*) 'vel_interpolation: cur_n_x_A = ', cur_n_x_A, ' cur_n_y_A = ', cur_n_y_A, ' cur_n_z_A = ', cur_n_z_A
 ! write(*,*) 'vel_interpolation: n_A = ', n_A

 IF(cur_x > xmin + w_vgrid_x / 2.0 .and. cur_x < xmax - w_vgrid_x /  2.0 .and.&
  & cur_y > ymin + w_vgrid_y / 2.0 .and. cur_y < ymax - w_vgrid_y /  2.0 .and. &
  & cur_z > zmin + w_vgrid_z / 2.0 .and. cur_z < zmax - w_vgrid_z /  2.0 ) THEN
  cur_n_x_B = floor((cur_pg_pos(ind_x) - xmin)/w_vgrid_x - 1.0/2.0) + 1
  cur_n_y_B = floor((cur_pg_pos(ind_y) - ymin)/w_vgrid_y - 1.0/2.0) + 1
  cur_n_z_B = floor((cur_pg_pos(ind_z) - zmin)/w_vgrid_z - 1.0/2.0) + 1
  n_B = cur_n_x_B + (N_vgrid_x - 1) * (cur_n_y_B - 1) + (N_vgrid_x - 1) * (N_vgrid_y - 1) * (cur_n_z_B - 1)
  cur_center_B = (/ w_vgrid_x * (cur_n_x_B + 1.D0) + xmin, w_vgrid_y * (cur_n_y_B + 1.D0) + ymin, &
                   & w_vgrid_z * (cur_n_z_B + 1.D0) + zmin /)
 ELSE
  n_B = 0
 END IF

 ! what is the best grid, A or B?
 dist_A = sqrt((cur_pg_pos(ind_x) - cur_center_A(ind_x))**2 + (cur_pg_pos(ind_y) - cur_center_A(ind_y))**2 + &
   & (cur_pg_pos(ind_z) - cur_center_A(ind_z))**2)
 ! write(*,*) 'vel_interpolation: dist_A = ', dist_A/R_star
 IF(n_b > 0) THEN
  dist_B = sqrt((cur_pg_pos(ind_x) - cur_center_B(ind_x))**2 + (cur_pg_pos(ind_y) - cur_center_B(ind_y))**2 + &
   &(cur_pg_pos(ind_y) - cur_center_B(ind_y))**2)
 ELSE IF (n_B == 0) THEN
  dist_B = large_number
 END IF
  

 !_______________________________________________________________
 !               CHOISE OF MODGRID POINTS
 !_______________________________________________________________
 ! choosing the correct modGrid points
 ! this is done for unification of the forthcoming code (after this if)
 IF(dist_A < dist_B) THEN
  cur_n_points = n_points_A(n_A)
  cur_start_index = indices_A(n_A)
  cur_end_index = cur_start_index + cur_n_points - 1

  ALLOCATE(cur_points(cur_n_points))

  cur_points = vg_indexy_A(cur_start_index:cur_end_index,2)
  chosen_gridAB = grid_A
  
 ELSE
  cur_n_points = n_points_B(n_B)
  cur_start_index = indices_B(n_B)
  cur_end_index = cur_start_index + cur_n_points - 1

  ! write(*,*) 'vel_interpolation: n_closest = ', n_closest
  ! interp_dist(I,J)
  ! interp_dist(:,1) -- distance
  ! interp_dist(:,2) -- index of the point
  ALLOCATE(cur_points(cur_n_points))

  cur_points = vg_indexy_B(cur_start_index:cur_end_index,2)
  DO cur_vg_point = 1, cur_n_points
   cur_vg_index = vg_indexy_B(cur_vg_point,1)
  END DO
  chosen_gridAB = grid_B
 END IF ! dist_A < dist_B

 ! DO cur_index_i = 1, SIZE(cur_points)
 !  write(*,*) 'vel_interpolation: pos = ', model_grid(cur_index_i)%vec_pos
 ! END DO
 ! write(*,*) 'vel_interpolation: cur_start_index = ', cur_start_index, ' cur_end_index = ', cur_end_index
 ! write(*,*) 'vel_interpolation: cur_n_points = ', cur_n_points

 
 !_______________________________________________________________
 !            SEEKING FOR THE n_closest CLOSEST POINTS
 !_______________________________________________________________
 interp_dist(:,ind_dist) = large_number
 interp_dist(:,ind_index) = 0.D0
 seeking = .TRUE.
 ! going through the array interp_dist and sorting from the smallest distances to largest distances
 DO cur_nearest_point = 1, cur_n_points ! #L00
  IF(procout) write(*,*) 'vel_interpolation: _______________________________________________________________________'
  IF(procout) write(*,*) 'vel_interpolation: cur_point = ', cur_nearest_point, '/', cur_n_points, ' dist = ', dist
  cur_nop = cur_points(cur_nearest_point)
  cur_mg_pos = model_grid(cur_nop)%vec_pos
  ! write(45,*) cur_mg_pos
  ! distance of a current modGrid point and the 
  dist = sqrt((cur_pg_pos(ind_x) - cur_mg_pos(ind_x))**2 + (cur_pg_pos(ind_y) - cur_mg_pos(ind_y))**2 + &
   & (cur_pg_pos(ind_z) - cur_mg_pos(ind_z))**2)
  write(*,*) 'vel_interpolation pg = ', cur_pg_pos/R_star, 'mg = ', cur_mg_pos/R_star

  ! the current distance in the array is larger than the current distance, we should move all values
  ! one level below and save the new distance to the current index
  !write(*,*) 'vel_interpolation: intpdist = ', interp_dist(cur_index_pos, ind_dist), ' dist = ', dist
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
  write(*,*) 'vel_interpolation: cur_mg_pos = ', cur_mg_pos

  ! if the current point is not closer than the most far (last index) point in the interp_dist,
  ! we do not have to do anything and go to the following point
  IF(interp_dist(n_closest, ind_dist) < dist) THEN
   IF(procout) write(*,*) 'vel_interpolation: the point ', cur_nearest_point, ' is too far away'
   CYCLE
  END IF

  ! in this case, we have to go through the saved list of points to see, whether we can save the new point or not
  ! IF(n_coor(coor_x) >= 2 .or. n_coor(coor_y) >= 2 .or. n_coor(coor_z) >= 2) THEN
  count_xyz(:) = 0
  IF(procout) write(*,*) 'vel_interpolation: seeking for the points on the same line'
  cur_del_index = 1
  DO cur_index_I = 1, n_closest
   cur_iti_mgi = INT(interp_dist(cur_index_I, ind_index))
   DO cur_index_J = cur_index_I + 1, n_closest
    cur_itj_mgi = INT(interp_dist(cur_index_J, ind_index))
    IF(INT(cur_iti_mgi) /= 0 .and. INT(cur_itj_mgi) /= 0) THEN
     IF(procout) write(*,*) 'vel_interpolation: index = ', interp_dist(cur_index_I, ind_index), interp_dist(cur_index_J, ind_index)
     cur_saved_mg_pos_A = model_grid(cur_iti_mgi)%vec_pos
     cur_saved_mg_pos_B = model_grid(cur_itj_mgi)%vec_pos

     vec_AB = (cur_saved_mg_pos_A - cur_saved_mg_pos_B)/NORM2(cur_saved_mg_pos_A - cur_saved_mg_pos_B)
     vec_AC = (cur_saved_mg_pos_A - cur_mg_pos)/NORM2(cur_saved_mg_pos_A - cur_mg_pos)
     IF(procout) write(*,*) 'vel_interpolation: vektory vec_AB a vec_AC'
     IF(procout) write(*,*) 'vec_AB = ', vec_AB
     IF(procout) write(*,*) 'vec_AC = ', vec_AC
     ! testing if a new point is on the same line as the other points already saved into interp_dist
     IF(((vec_AB(ind_x) == vec_AC(ind_x)) .and. (vec_AB(ind_y) == vec_AC(ind_y)) .and. (vec_AB(ind_z) == vec_AC(ind_z))) .or. &
      & ((vec_AB(ind_x) == -vec_AC(ind_x)) .and. (vec_AB(ind_y) == -vec_AC(ind_y)) .and. (vec_AB(ind_z) == -vec_AC(ind_z)))) THEN 
      ! zatím to risknu a vezmu ten druhý index, který by měl být index vzdálenějšího bodu
      f_indexy(cur_del_index) = cur_index_J
      mgi_indexy(cur_del_index) = cur_itj_mgi
      cur_del_index = cur_del_index + 1
      IF(procout) write(*,*) 'vel_interpolation: for delete = ', cur_index_J
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
    IF(f_indexy(1) > 0) index_delete(ind_x) = f_indexy(ind_x)
    IF(f_indexy(2) > 0) index_delete(ind_y) = f_indexy(ind_y)
    IF(f_indexy(3) > 0) index_delete(ind_z) = f_indexy(ind_z)
    IF(procout) write(*,*) 'vel_interpolation: index_delete = ', index_delete(:)
    nahrada = .true.
   ELSE
    nahrada = .false.
   END IF
   IF(procout) write(*,*) 'vel_interpolation: index to delete: ', index_delete(:)

   ! we now delete the old point and add a better point into the array
   ! deleting a point on a position index_delete
   DO cur_index_K = 1, const_dimofspace
    IF(index_delete(cur_index_K) > 0) THEN
     if(procout) write(*,*) 'vel_interpolation: smazani bodu ', index_delete(cur_index_K), ' mod_grid = ',&
      & interp_dist(index_delete(cur_index_K), ind_index)
     DO cur_index_I = 1, n_closest
      cur_mgi = INT(interp_dist(cur_index_I, ind_index))
      IF(cur_mgi == mgi_indexy(cur_index_K)) THEN
       IF(procout) write(*,*) 'vel_interpolation: smazani indexu: ', cur_mgi
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
  write(*,*) 'vel_interpolation: novyBod = ', novyBod
  IF(novyBod .or. nahrada) THEN
   pom = interp_dist
   write(*,*) 'vel_interpolation: cur_index_pos = ', cur_index_pos, ' n_closest = ', n_closest
   DO cur_index_pos = 1, n_closest
    write(*,*) 'vel_interpolation: cip = ', cip, ' cur_index_pos = ', cur_index_pos
    cur_dist = interp_dist(cur_index_pos, ind_dist)
    write(*,*) 'vel_interpolation: cur_dist = ', cur_dist, ' dist = ', dist
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
  write(*,*) 'vel_interpolation: interp_dist() = ', INT(interp_dist(:,ind_index))
  !IF(.NOT. seeking) EXIT
  ! CALL vel_intp_choice(interp_dist, n_closest, 2, cur_pg_pos)
 END DO ! #L00 loop over all possible points


 DO cur_index_J = 1, n_closest
  cur_iti_mgi = INT(interp_dist(cur_index_J, ind_index))
  IF(cur_iti_mgi == 0) THEN
   DO cur_nearest_point = 1, cur_n_points
    cur_nop = cur_points(cur_nearest_point)
    cur_mg_pos = model_grid(cur_nop)%vec_pos
    ! write(41,*) cur_mg_pos
   END DO
   EXIT
  ELSE
   ! write(42,*) model_grid(cur_iti_mgi)%vec_pos
  END IF
 END DO
 ! ARTIFICIAL POINTS
 ! if the interp_dist still contains zeros, we have to add some artificial points
 ! TOTO JEŠTĚ NENÍ DOKONČENO!!!
 DO cur_index_J = 1, n_closest
  cur_iti_mgi = INT(interp_dist(cur_index_J, ind_index))
  IF(cur_iti_mgi == 0) THEN

  END IF
 END DO

 ! write(*,*) 'vel_interpolation: interp_dist(:,2) = ', interp_dist(:,2)
 ! write(40,*) cur_pg_pos
 CALL vel_vector_interpolation(cur_pg_pos, interp_dist(:,2), n_closest, vel_vector)
 ! write(*,*) 'vel_interpolation:************************************************'
 ! write(*,*) 'vel_interpolation: interp_dist = ', interp_dist(:,ind_index)

 IF(isnan(vel_vector(ind_x)) .or. isnan(vel_vector(ind_y)) .or. isnan(vel_vector(ind_z))) THEN 
  DO cur_index_I = 1, cur_n_points
   cur_index = cur_points(cur_index_I)
   cur_mg_pos = model_grid(cur_index)%vec_pos
   ! write(45,*) cur_mg_pos
  END DO
  DO  cur_index_I = 1, n_closest
   cur_mgi = INT(interp_dist(cur_index_I, ind_index))
   cur_mg_pos = model_grid(cur_mgi)%vec_pos
   write(43,*) cur_mg_pos, vel_vector, INT(interp_dist(cur_index_I, ind_index))
  END DO
  write(44,*) cur_pg_pos
  STOP 'vel_vector_interpolation: NaNs'
 END IF

 DEALLOCATE(cur_points)

 dyn_cell(cur_prop_cell)%vec_vel = vel_vector
 ! write(43,*) cur_pg_pos, vel_vector

END DO ! a loop over all propGrid cells







END SUBROUTINE vel_interpolation
