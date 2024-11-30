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
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_center_A, cur_center_B

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
CALL get_scalar_index(cur_pos, .true., n_A)
CALL get_scalar_index(cur_pos, .false., n_B)
                 ! & w_vgrid_z * (cur_n_z_A + 5.D-1) + zmin /)

! write(*,*) 'vel_interpolation: cur_prop_cell = ', cur_prop_cell
!_______________________________________________________________
!  #              CHOISE OF GRID A OR B
!_______________________________________________________________
! the index in the AB grid
CALL get_vg_centre_A(cur_pos, cur_center_A)
CALL get_vg_centre_B(cur_pos, cur_center_B)

! what is the best grid, A or B?
dist_A = sqrt((cur_pos(ind_x) - cur_center_A(ind_x))**2 + (cur_pos(ind_y) - cur_center_A(ind_y))**2)
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
 ! write(*,*) 'n_closest_points_2D: n_A = ', n_A
 cur_n_points = n_points_A(n_A)
 ALLOCATE(cur_points(cur_n_points))
 IF(cur_n_points > 0) THEN
  cur_start_index = indices_A(n_A)
  cur_end_index = cur_start_index + cur_n_points - 1
  cur_points = vg_indexy_A(cur_start_index:cur_end_index,ind_mg)
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
 cur_points = vg_indexy_B(cur_start_index:cur_end_index,ind_mg)
 DO cur_vg_point = 1, cur_n_points
  cur_vg_index = vg_indexy_B(cur_vg_point,ind_vg)
 END DO
 chosen_gridAB = grid_B
END IF ! dist_A < dist_B

CALL seek_nclosest_points(cur_pos, cur_points, N_points, cur_n_points, closest_points)

! write(*,*) 'n_closest_points_2D: cur_points = ', cur_points

! write(*,*) 'n_closest_points_2D: closest_points = ', closest_points


END SUBROUTINE n_closest_points_2D
