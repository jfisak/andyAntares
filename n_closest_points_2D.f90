! the subroutine for obtaining N_points modGrid points (cells) closest to a particular
! modGrid point
! 
! INPUT: N_points, INT -- number of points
! OUTPUT: closest_points, INT(N_points) -- returns indeces of the points
!
SUBROUTINE n_closest_points_2D(cur_mgi, N_points, closest_points)


USE types
USE virt_gridAB
IMPLICIT NONE

INTEGER                                 :: N_points
INTEGER, DIMENSION(N_points)            :: closest_points
INTEGER, PARAMETER                      :: min_incell = 10

DOUBLE PRECISION, PARAMETER             :: large_number = 1.D90

INTEGER, ALLOCATABLE                    :: cur_points(:)
INTEGER                                 :: cur_start_index, cur_end_index, cur_n_points
DOUBLE PRECISION                        :: dist_A, dist_B! , dist_C

INTEGER, PARAMETER                      :: n_closest = 8

INTEGER                                 :: chosen_gridAB
INTEGER, PARAMETER                      :: grid_A = 1, grid_B = 2
INTEGER, PARAMETER                      :: ind_dist = 1, ind_index = 2
DOUBLE PRECISION, DIMENSION(2)           :: cur_center_A, cur_center_B

INTEGER                                 :: cur_vg_index, cur_vg_point

DOUBLE PRECISION                        :: cur_x, cur_y, cur_z
INTEGER                                 :: cur_n_x_A, cur_n_y_A!, cur_n_z_A
DOUBLE PRECISION                        :: cur_n_x_B, cur_n_y_B! , cur_n_z_B
INTEGER                                 :: n_A, n_B

INTEGER, PARAMETER                                      :: coor_x = 1, coor_y = 2, coor_z = 3

INTEGER                                                 :: cur_mgi

DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_pos

cur_x = model_grid(cur_mgi)%rwind
cur_y = model_grid(cur_mgi)%angle
cur_z = 0.D0

cur_pos = (/ cur_x, cur_y, cur_z/)
! calculation of the current grid indeces
cur_n_x_A = floor((cur_pos(ind_x) - vg_xmin)/w_vgrid_x) + 1
cur_n_y_A = floor((cur_pos(ind_y) - vg_ymin)/w_vgrid_y) + 1
n_A = INT(cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1))! + N_vgrid_x * N_vgrid_y * (cur_n_z_A - 1)
cur_center_A = (/ w_vgrid_x * (cur_n_x_A + 5.D-1) + vg_xmin, w_vgrid_y * (cur_n_y_A + 5.D-1) + vg_ymin /)!, &
                 ! & w_vgrid_z * (cur_n_z_A + 5.D-1) + zmin /)

! write(*,*) 'vel_interpolation: cur_prop_cell = ', cur_prop_cell
!_______________________________________________________________
!  #              CHOISE OF GRID A OR B
!_______________________________________________________________
! the index in the AB grid
n_A = INT(cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1))! + N_vgrid_x * N_vgrid_y * (cur_n_z_A - 1)
cur_center_A = (/ w_vgrid_x * (cur_n_x_A + 5.D-1) + vg_xmin, w_vgrid_y * (cur_n_y_A + 5.D-1)/)
! write(*,*) 'n_closest_points_2D: cur_pos = ', cur_pos
! write(*,*) 'n_closest_points_2D: w_vgrid_x = ', w_vgrid_x, ' w_vgrid_y = ', w_vgrid_y
! write(*,*) 'n_closest_points_2D: vg_xmin = ', vg_xmin, ' vg_ymin = ', vg_ymin
! write(*,*) 'n_closest_points_2D: cur_n_x_A = ', cur_n_x_A, ' cur_n_y_A = ', cur_n_y_A
! write(*,*) 'n_closest_points_2D: n_A = ', n_A

IF(cur_pos(ind_x) > vg_xmin + w_vgrid_x / 2.0 .and. cur_pos(ind_x) < vg_xmax - w_vgrid_x /  2.0 .and.&
 & cur_pos(ind_y) > vg_ymin + w_vgrid_y / 2.0 .and. cur_pos(ind_y) < vg_ymax - w_vgrid_y /  2.0) THEN
 cur_n_x_B = floor((cur_pos(ind_x) - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
 cur_n_y_B = floor((cur_pos(ind_y) - vg_ymin)/w_vgrid_y - 1.0/2.0) + 1
 n_B = INT(cur_n_x_B + (N_vgrid_x - 1) * (cur_n_y_B - 1))! + (N_vgrid_x - 1) * (N_vgrid_y - 1) * (cur_n_z_B - 1)
 cur_center_B = (/ w_vgrid_x * (cur_n_x_A + 1.D0) + vg_xmin, w_vgrid_y * (cur_n_y_A + 1.D0) + vg_ymin /)!, &
                  ! & w_vgrid_z * (cur_n_z_A + 1.D0) + vg_zmin /)
ELSE
 n_B = 0
END IF

! what is the best grid, A or B?
dist_A = sqrt((cur_pos(ind_x) - cur_center_A(1))**2 + (cur_pos(ind_y) - cur_center_A(2))**2)
! write(*,*) 'vel_interpolation: dist_A = ', dist_A/R_star
IF(n_b > 0) THEN
 dist_B = sqrt((cur_pos(ind_x) - cur_center_B(1))**2 + (cur_pos(ind_y) - cur_center_B(2))**2)
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
 ALLOCATE(cur_points(cur_n_points))
 IF(cur_n_points > 0) THEN
  cur_start_index = indices_A(n_A)
  cur_end_index = cur_start_index + cur_n_points - 1
  cur_points = vg_indexy_A(cur_start_index:cur_end_index,2)
  chosen_gridAB = grid_A
 ELSE
  cur_start_index = -1
  cur_end_index = -1
 END IF

 
ELSE
 cur_n_points = n_points_B(n_B)
 ALLOCATE(cur_points(cur_n_points))
 cur_start_index = indices_B(n_B)
 cur_end_index = cur_start_index + cur_n_points - 1

 ! write(*,*) 'vel_interpolation: n_closest = ', n_closest

 ! write(*,*) 'n_closest_points_2D: vg_indexy_B = ', vg_indexy_B(:,:)
 cur_points = vg_indexy_B(cur_start_index:cur_end_index,2)
 DO cur_vg_point = 1, cur_n_points
  cur_vg_index = vg_indexy_B(cur_vg_point,1)
 END DO
 chosen_gridAB = grid_B
END IF ! dist_A < dist_B

CALL seek_nclosest_points(cur_pos, cur_points, N_points, cur_n_points, closest_points)

! write(*,*) 'n_closest_points_2D: cur_points = ', cur_points

! write(*,*) 'n_closest_points_2D: closest_points = ', closest_points


END SUBROUTINE n_closest_points_2D
