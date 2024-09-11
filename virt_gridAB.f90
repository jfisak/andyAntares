MODULE virt_gridAB

USE types

IMPLICIT NONE


INTEGER                                 :: N_vgrid_x, N_vgrid_y, N_vgrid_z
INTEGER                                 :: N_vgrid_cells_A, N_vgrid_cells_B
INTEGER, ALLOCATABLE                    :: n_points_A(:), n_points_B(:), &
                                         & indices_A(:), indices_B(:)
DOUBLE PRECISION                        :: w_vgrid_x, w_vgrid_y, w_vgrid_z

INTEGER, ALLOCATABLE                    :: vg_indexy_A(:,:), vg_indexy_B(:,:)


DOUBLE PRECISION                        :: vg_xmin, vg_xmax, vg_ymin, vg_ymax, vg_zmin, vg_zmax

LOGICAL                                 :: is_initialized = .false.

CONTAINS

SUBROUTINE virt_gridAB_init()

 USE types
 USE constants
 IMPLICIT NONE

INTEGER, DIMENSION(2)                   :: dummy_var_A, dummy_var_B, dummy_A, dummy_B

INTEGER                                 :: cur_ind_A, cur_ind_B, cur_vpg_cell
INTEGER                                 :: cur_point
INTEGER                                 :: n_A, n_B
INTEGER                                 :: n_in_cell
INTEGER                                 :: n_zeros
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_center_A, cur_center_B
INTEGER, PARAMETER                      :: min_incell = 10

INTEGER                                 :: cur_iter

DOUBLE PRECISION                        :: cur_x, cur_y, cur_z
INTEGER                                 :: cur_n_x_A, cur_n_y_A, cur_n_z_A
INTEGER                                 :: cur_n_x_B, cur_n_y_B, cur_n_z_B

INTEGER                                 :: cur_vgi, cur_ind_sorted, cur_mgi

INTEGER, ALLOCATABLE                    :: counter_A(:), counter_B(:)

 n_in_cell = FLOOR((n_modelgrid/min_incell)**(1.0/3.0))
 IF(n_in_cell < 1) n_in_cell = 1
 
 N_vgrid_x = n_in_cell
 N_vgrid_y = n_in_cell
 N_vgrid_z = n_in_cell
 N_vgrid_cells_A = N_vgrid_x * N_vgrid_y * N_vgrid_z
 N_vgrid_cells_B = (N_vgrid_x - 1) * (N_vgrid_y - 1) * (N_vgrid_z - 1)
 ALLOCATE(n_points_A(N_vgrid_cells_A), n_points_B(N_vgrid_cells_B))
 ALLOCATE(indices_A(N_vgrid_cells_A), indices_B(N_vgrid_cells_B))
 ALLOCATE(vg_indexy_A(n_modelgrid, 2), vg_indexy_B(n_modelgrid, 2))
 n_points_A(:) = 0
 n_points_B(:) = 0
 n_zeros = 0
 
 IF(model_type == 1) THEN
  vg_xmin = R_star
  vg_xmax = R_inf
  w_vgrid_x = (R_inf - R_star)/N_vgrid_x
  vg_ymin = ymin
  vg_ymax = ymax
  w_vgrid_y = 0.D0
  vg_zmin = zmin
  vg_zmax = zmax
  w_vgrid_z = 0.D0
 END IF
 IF(model_type == 2) THEN
  IF(inputmodel == 1) THEN
   vg_xmin = R_star
   vg_xmax = R_inf
   w_vgrid_x = (R_inf - R_star)/N_vgrid_x
   vg_ymin = 0.D0
   vg_ymax = 2.D0 * const_pi
   w_vgrid_y = const_pi/N_vgrid_y
   vg_zmin = zmin
   vg_zmax = zmax
   w_vgrid_z = 0.D0
  END IF
 END IF
 IF(model_type == 3) THEN
  vg_xmin = xmin
  vg_xmax = xmax
  w_vgrid_x = abs(vg_xmax - vg_xmin)/N_vgrid_x
  vg_ymin = ymin
  vg_ymax = ymax
  w_vgrid_y = abs(vg_ymax - vg_ymin)/N_vgrid_y
  vg_zmin = zmin
  vg_zmax = zmax
  w_vgrid_z = abs(vg_zmax - vg_zmin)/N_vgrid_z
 END IF

 !_______________________________________________________________
 !   #01        CALCULATE VG INDEX OF MG CELLS
 !_______________________________________________________________
 ! calculation of the virGrid index
 DO cur_point = 1, n_modelgrid
  IF(model_grid(cur_point)%assoc_cells == 0) CYCLE

  IF(model_type == 1) THEN
   cur_x = model_grid(cur_point)%rwind
   cur_y = 0.D0
   cur_z = 0.D0
  ELSE IF(model_type == 2) THEN
   IF(inputmodel == 1) THEN
    cur_x = model_grid(cur_point)%rwind
    cur_y = model_grid(cur_point)%angle
    cur_z = 0.D0
   END IF
  END IF
  IF(model_type == 3) THEN
   cur_x = model_grid(cur_point)%vec_pos(ind_x)
   cur_y = model_grid(cur_point)%vec_pos(ind_y)
   cur_z = model_grid(cur_point)%vec_pos(ind_z)
  END IF

  IF(model_type == 1) THEN
   cur_n_x_A = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
   cur_n_y_A = 1
   cur_n_z_A = 1
  ELSE IF(model_type == 2) THEN
   cur_n_x_A = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
   cur_n_y_A = floor((cur_y - vg_ymin)/w_vgrid_y) + 1
   cur_n_z_A = 1
  ELSE IF(model_type == 3) THEN
   cur_n_x_A = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
   cur_n_y_A = floor((cur_y - vg_ymin)/w_vgrid_y) + 1
   cur_n_z_A = floor((cur_z - vg_zmin)/w_vgrid_z) + 1
  END IF
  
  cur_center_A = (/ w_vgrid_x * (cur_n_x_A + 5.D-1) + vg_xmin, w_vgrid_y * (cur_n_y_A + 5.D-1) + vg_ymin, &
                  & w_vgrid_z * (cur_n_z_A + 5.D-1) + vg_zmin /)
  ! write(*,*) 'vel_interpolation: cur_center_A = ', cur_center_A
 
  ! n_A -- numerical index of VG cell
  IF(model_type == 1) THEN
   n_A = cur_n_x_A
  ELSE IF(model_type == 2) THEN
   n_A = cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1)
  ELSE IF(model_type == 3) THEN
   n_A = cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1) + N_vgrid_x * N_vgrid_y * (cur_n_z_A - 1)
  END IF
  ! write(*,*) 'vel_interpolation: n_x = ', cur_n_x_A, ' n_y = ', cur_n_y_A, ' n_z = ', cur_n_z_A
  ! write(*,*) 'vel_interpolation: n_A = ', n_A
  ! n_points -- number of points for the given cell
  n_points_A(n_A) = n_points_A(n_A) + 1
  ! vg_indexy -- list of indeces model grid --> VG index point
  vg_indexy_A(cur_point, 1) = n_A
  vg_indexy_A(cur_point, 2) = cur_point
  !!!!!!!!
  ! repete for the B grid
  IF(cur_x > vg_xmin + w_vgrid_x/2.0 .and. cur_x < vg_xmax - w_vgrid_x/2.0 .and.&
   & cur_y > vg_ymin + w_vgrid_y/2.0 .and. cur_y < vg_ymax - w_vgrid_y/2.0 .and. &
   & cur_z > vg_zmin + w_vgrid_z/2.0 .and. cur_z < vg_zmax - w_vgrid_z/2.0 ) THEN
   IF(model_type == 1) THEN
    cur_n_x_B = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
    cur_n_y_B = 1
    cur_n_z_B = 1
   ELSE IF(model_type == 2) THEN
    cur_n_x_B = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
    cur_n_y_B = floor((cur_y - vg_ymin)/w_vgrid_y - 1.0/2.0) + 1
    cur_n_z_B = 1
   ELSE IF(model_type == 3) THEN
    cur_n_x_B = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
    cur_n_y_B = floor((cur_y - vg_ymin)/w_vgrid_y - 1.0/2.0) + 1
    cur_n_z_B = floor((cur_z - vg_zmin)/w_vgrid_z - 1.0/2.0) + 1
   END IF
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
 write(*,*) 'virt_gridAB_init: started sorting A'
 ! calculating number of points in each vg grid cell (A)
 ! now we know where the indeces will be located, we set up an array
 ! containing the initial indeces for all vg indeces
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

 counter_A = n_points_A
 counter_B = n_points_B

 DO cur_mgi = 1, n_modelgrid
  ! A
  cur_vgi = vg_indexy_A(cur_mgi, 1)
  IF(counter_A(cur_vgi) > 0) THEN
   cur_ind_sorted = indices_A(cur_vgi) + counter_A(cur_vgi) - 1
  ELSE
   write(*,*) 'virt_gridAB_init: sorting is not OK'
   write(*,*) 'want to add a point of cur_vgi = ', cur_vgi
   write(*,*) 'with no point left'
   STOP 'virt_gridAB_init'
  END IF
  counter_A(cur_vgi) = counter_A(cur_vgi) - 1
  ! B
  cur_vgi = vg_indexy_B(cur_mgi, 1)
  IF(counter_B(cur_vgi) > 0) THEN
   cur_ind_sorted = indices_B(cur_vgi) + counter_B(cur_vgi) - 1
  ELSE
   write(*,*) 'virt_gridAB_init: sorting is not OK'
   write(*,*) 'want to add a point of cur_vgi = ', cur_vgi
   write(*,*) 'with no point left'
   STOP 'virt_gridAB_init'
  END IF
  counter_B(cur_vgi) = counter_B(cur_vgi) - 1
 END DO


 write(*,*) 'virt_gridAB_init: ended sorting A'
 
 !_______________________________________________________________
 ! sort the vg_indexy according to the VG index
 write(*,*) 'virt_gridAB_init: started sorting B'
 !_______________________________________________________________
 write(*,*) 'virt_gridAB_init: ended sorting B'
 
 !_______________________________________________________________
 ! index array
 !_______________________________________________________________
 

 is_initialized = .true.

END SUBROUTINE virt_gridAB_init

END MODULE virt_gridAB
