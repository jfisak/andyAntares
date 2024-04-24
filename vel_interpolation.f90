! this sbr will interpolate the vector field of velocity from the modGrid onto the propGrid
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

DOUBLE PRECISION, DIMENSION(3)          :: cur_pg_pos, cur_mg_pos
INTEGER                                 :: cur_prop_cell

DOUBLE PRECISION, PARAMETER             :: large_number = 1.D90

INTEGER, ALLOCATABLE                    :: cur_points(:)
DOUBLE PRECISION, DIMENSION(3)          :: cur_center_A, cur_center_B
INTEGER                                 :: cur_start_index, cur_end_index, cur_n_points
DOUBLE PRECISION                        :: dist_A, dist_B

INTEGER, PARAMETER                      :: n_closest = 8
DOUBLE PRECISION, ALLOCATABLE           :: interp_dist(:,:), pom(:,:)

INTEGER                                 :: chosen_gridAB
INTEGER, PARAMETER                      :: grid_A = 1, grid_B = 2
INTEGER, PARAMETER                      :: ind_dist = 1, ind_index = 2

INTEGER                                 :: cur_nearest_point, cur_nop, cur_vg_index, cur_vg_point
INTEGER                                 :: cip, cur_index, cur_index_pos
DOUBLE PRECISION                        :: dist
LOGICAL                                 :: seeking
! INTEGER                                 :: cur_index_i

DOUBLE PRECISION, DIMENSION(3)          :: vel_vector

INTEGER                                 :: n_coor_x, n_coor_y, n_coor_z
INTEGER                                 :: count_x, count_y, count_z
DOUBLE PRECISION, DIMENSION(3)          :: cur_saved_mg_pos
LOGICAL                                 :: novyBod
INTEGER                                 :: cur_iti_mgi
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
 cur_x = model_grid(cur_point)%vec_pos(1)
 cur_y = model_grid(cur_point)%vec_pos(2)
 cur_z = model_grid(cur_point)%vec_pos(3)
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
 cur_n_x_A = floor((cur_pg_pos(1) - xmin)/w_vgrid_x) + 1
 cur_n_y_A = floor((cur_pg_pos(2) - ymin)/w_vgrid_y) + 1
 cur_n_z_A = floor((cur_pg_pos(3) - zmin)/w_vgrid_z) + 1
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
  cur_n_x_B = floor((cur_pg_pos(1) - xmin)/w_vgrid_x - 1.0/2.0) + 1
  cur_n_y_B = floor((cur_pg_pos(2) - ymin)/w_vgrid_y - 1.0/2.0) + 1
  cur_n_z_B = floor((cur_pg_pos(3) - zmin)/w_vgrid_z - 1.0/2.0) + 1
  n_B = cur_n_x_B + (N_vgrid_x - 1) * (cur_n_y_B - 1) + (N_vgrid_x - 1) * (N_vgrid_y - 1) * (cur_n_z_B - 1)
  cur_center_B = (/ w_vgrid_x * (cur_n_x_A + 1.D0) + xmin, w_vgrid_y * (cur_n_y_A + 1.D0) + ymin, &
                   & w_vgrid_z * (cur_n_z_A + 1.D0) + zmin /)
 ELSE
  n_B = 0
 END IF

 ! what is the best grid, A or B?
 dist_A = sqrt((cur_pg_pos(1) - cur_center_A(1))**2 + (cur_pg_pos(2) - cur_center_A(2))**2 + &
   & (cur_pg_pos(3) - cur_center_A(3))**2)
 ! write(*,*) 'vel_interpolation: dist_A = ', dist_A/R_star
 IF(n_b > 0) THEN
  dist_B = sqrt((cur_pg_pos(1) - cur_center_B(1))**2 + (cur_pg_pos(2) - cur_center_B(2))**2 + &
   &(cur_pg_pos(3) - cur_center_B(3))**2)
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
 ! write(*,*) 'vel_interpolation: interp_dist = ', interp_dist(:,ind_dist)
 DO cur_nearest_point = 1, cur_n_points
  cur_nop = cur_points(cur_nearest_point)
  cur_mg_pos = model_grid(cur_nop)%vec_pos
  ! write(45,*) cur_mg_pos
  ! write(*,*) 'vel_interpolation: cur_index = ', cur_nearest_point, ' cur_mg_pos = ', cur_mg_pos
  
  ! distance of a current modGrid point and the 
  dist = sqrt((cur_pg_pos(1) - cur_mg_pos(1))**2 + (cur_pg_pos(2) - cur_mg_pos(2))**2 + &
   & (cur_pg_pos(3) - cur_mg_pos(3))**2)
  ! write(*,*) 'vel_interpolation pg = ', cur_pg_pos/R_star, 'mg = ', cur_mg_pos/R_star
  ! write(*,*) 'vel_interpolation: cur_point = ', cur_nearest_point, ' dist = ', dist

  ! going through the array interp_dist and sorting from the smallest distances to largest distances
  DO cur_index_pos = 1, n_closest
   ! the current distance in the array is larger than the current distance, we should move all values
   ! one level below and save the new distance to the current index
   !write(*,*) 'vel_interpolation: intpdist = ', interp_dist(cur_index_pos, ind_dist), ' dist = ', dist
   seeking = .TRUE.
   IF(interp_dist(cur_index_pos, ind_dist) > dist) THEN
    ! testing of the coordinates
    n_coor_x = 0
    n_coor_y = 0
    n_coor_z = 0
    count_x = 0
    count_y = 0
    count_z = 0
    novyBod = .true.
    ! going value by value and counting the number of points with very same coordinate
    ! the number larger than two is useless, therefore we skip this point
    DO cip = 1, n_closest
     IF(interp_dist(cip, ind_dist) < large_number .and. INT(interp_dist(cip, ind_index)) /= 0) THEN
      cur_iti_mgi = INT(interp_dist(cip, ind_index))
      cur_saved_mg_pos = model_grid(cur_iti_mgi)%vec_pos
      write(*,*) 'vel_interpolation: cur_saved_mg_pos = ', cur_saved_mg_pos
      IF(cur_mg_pos(1) == cur_saved_mg_pos(1)) n_coor_x = n_coor_x + 1
      IF(cur_mg_pos(2) == cur_saved_mg_pos(2)) n_coor_y = n_coor_y + 1
      IF(cur_mg_pos(3) == cur_saved_mg_pos(3)) n_coor_z = n_coor_z + 1
      write(*,*) 'vel_interpolation: n_x = ', n_coor_x, ' n_y = ', n_coor_y, ' n_z = ', n_coor_z
     END IF
    END DO

    ! x direction
    IF(n_coor_x >= 2) THEN
     write(*,*) 'vel_interpolation: n_coor_x = ', n_coor_x
     DO cip = 1, n_closest
      IF(INT(interp_dist(cip, ind_index)) /= 0) THEN
       write(*,*) 'vel_interpolation: index = ', interp_dist(cip, ind_index)
       cur_iti_mgi = INT(interp_dist(cip, ind_index))
       cur_saved_mg_pos = model_grid(cur_iti_mgi)%vec_pos
       IF(cur_mg_pos(1) == cur_saved_mg_pos(1)) THEN
        ! this is the second occurence of the same coordinate
        ! we can delete this point
        write(*,*) 'vel_interpolation: count_x = ', count_x
        IF(interp_dist(cip, ind_dist) > dist .and. count_x == 1) THEN
         ! we should test the second, more distant, point
         write(*,*) 'vel_interpolation: deleting a point'
         ! move the indeces one position up
         DO cur_index = cip, n_closest - 1
          interp_dist(cur_index,:) = interp_dist(cur_index + 1,:)
         END DO
         interp_dist(n_closest,ind_dist) = large_number
         interp_dist(n_closest, ind_index) = 0
        ELSE
         novyBod = .false.
        END IF
        count_x = count_x + 1
       END IF
      END IF ! interp_dist(mod_index) /= 0
     END DO
     n_coor_x = 2
     write(*,*) 'vel_interpolation: n_coor_x = ', n_coor_x
    END IF ! n_coor_x > 2

    ! y direction
    IF(n_coor_y >= 2) THEN
     write(*,*) 'vel_interpolation: n_coor_y = ', n_coor_y
     DO cip = 1, n_closest
      IF(INT(interp_dist(cip, ind_index)) /= 0) THEN
       cur_iti_mgi = INT(interp_dist(cip, ind_index))
       cur_saved_mg_pos = model_grid(cur_iti_mgi)%vec_pos
       IF(cur_mg_pos(2) == cur_saved_mg_pos(2)) THEN
        ! this is the second occurence of the same coordinate
        ! we can delete this point
        write(*,*) 'vel_interpolation: count_y = ', count_y
        IF(interp_dist(cip, ind_dist) > dist .and. count_y == 1) THEN
         ! we should test the second, more distant, point
         write(*,*) 'vel_interpolation: deleting a point'
         ! move the indeces one position up
         DO cur_index = cip, n_closest - 1
          interp_dist(cur_index,:) = interp_dist(cur_index + 1,:)
         END DO
         interp_dist(n_closest,ind_dist) = large_number
         interp_dist(n_closest, ind_index) = 0
        ELSE
         novyBod = .false.
        END IF
        count_y = count_y + 1
       END IF
      END IF ! interp_dist(mod_index) /= 0
     END DO
     n_coor_y = 2
     write(*,*) 'vel_interpolation: n_coor_y = ', n_coor_y
    END IF ! n_coor_y > 2

    ! z direction
    IF(n_coor_z >= 2) THEN
     write(*,*) 'vel_interpolation: n_coor_z = ', n_coor_z
     DO cip = 1, n_closest
      IF(INT(interp_dist(cip, ind_index)) /= 0) THEN
       cur_iti_mgi = INT(interp_dist(cip, ind_index))
       cur_saved_mg_pos = model_grid(cur_iti_mgi)%vec_pos
       IF(cur_mg_pos(3) == cur_saved_mg_pos(3)) THEN
        ! this is the second occurence of the same coordinate
        ! we can delete this point
        write(*,*) 'vel_interpolation: count_z = ', count_z
        IF(interp_dist(cip, ind_dist) > dist .and. count_z == 1) THEN
         ! we should test the second, more distant, point
         write(*,*) 'vel_interpolation: deleting a point'
         ! move the indeces one position up
         DO cur_index = cip, n_closest - 1
          interp_dist(cur_index,:) = interp_dist(cur_index + 1,:)
         END DO
         interp_dist(n_closest,ind_dist) = large_number
         interp_dist(n_closest, ind_index) = 0
        ELSE
         novyBod = .false.
        END IF
        count_z = count_z + 1
       END IF
      END IF ! interp_dist(mod_index) /= 0
     END DO
     n_coor_z = 2
     write(*,*) 'vel_interpolation: n_coor_z = ', n_coor_z
    END IF ! n_coor_z > 2

    ! we can now put the point into the list
    write(*,*) 'vel_interpolation: novyBod = ', novyBod

    IF(novyBod) THEN
     pom = interp_dist
     DO cip = cur_index_pos, n_closest
      ! write(*,*) 'vel_interpolation: cip = ', cip, ' cur_index_pos = ', cur_index_pos
      IF(cip == cur_index_pos) THEN
       interp_dist(cip, ind_dist) = dist
       interp_dist(cip, ind_index) = cur_nop
      ELSE
       interp_dist(cip, ind_dist) = pom(cip - 1, ind_dist)
       interp_dist(cip, ind_index) = pom(cip - 1, ind_index)
      END IF
      seeking = .FALSE.
     END DO
    END IF ! novyBod
   END IF ! cur_dist < dist in array
   IF(.NOT. seeking) EXIT
  END DO ! cur_index_pos
  ! write(*,*) 'vel_interpolation: interp_dist = ', interp_dist(:,ind_dist)
 END DO
 
 ! DO cur_J = 1, n_closest
 !  IF(interp_dist(cur_J, ind_index) == 0) THEN
 !  END IF
 ! END DO
 ! write(*,*) 'vel_interpolation: interp_dist(:,2) = ', interp_dist(:,2)
 CALL vel_vector_interpolation(cur_pg_pos, interp_dist(:,2), n_closest, vel_vector)
 ! write(*,*) 'vel_interpolation:************************************************'
 ! write(*,*) 'vel_interpolation: interp_dist = ', interp_dist(:,ind_index)
 DEALLOCATE(cur_points)

 dyn_cell(cur_prop_cell)%vec_vel = vel_vector
 ! write(43,*) cur_pg_pos, vel_vector

END DO ! a loop over all propGrid cells







END SUBROUTINE vel_interpolation
