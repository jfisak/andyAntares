MODULE virt_gridAB

USE types

IMPLICIT NONE


INTEGER                                 :: N_vgrid_x, N_vgrid_y, N_vgrid_z
INTEGER                                 :: N_vgrid_cells_A, N_vgrid_cells_B
INTEGER, ALLOCATABLE                    :: n_points_A(:), n_points_B(:), &
                                         & indices_A(:), indices_B(:)
DOUBLE PRECISION                        :: w_vgrid_x, w_vgrid_y, w_vgrid_z

INTEGER, ALLOCATABLE                    :: vg_indexy_A(:,:), vg_indexy_B(:,:)
INTEGER, ALLOCATABLE                    :: vg_pom_A(:,:), vg_pom_B(:,:)


DOUBLE PRECISION                        :: vg_xmin, vg_xmax, vg_ymin, vg_ymax, vg_zmin, vg_zmax

LOGICAL                                 :: is_initialized = .false.

CONTAINS

SUBROUTINE virt_gridAB_init()

 USE types
 USE constants
 IMPLICIT NONE

INTEGER                                 :: cur_ind_A, cur_ind_B, cur_vpg_cell
INTEGER                                 :: cur_point
INTEGER                                 :: n_A, n_B
INTEGER                                 :: n_in_cell
INTEGER                                 :: n_zeros
INTEGER, PARAMETER                      :: min_incell = 10

INTEGER                                 :: cur_index_mgi

DOUBLE PRECISION                        :: cur_x, cur_y, cur_z
INTEGER                                 :: cur_n_x_A, cur_n_y_A, cur_n_z_A
INTEGER                                 :: cur_n_x_B, cur_n_y_B, cur_n_z_B

INTEGER                                 :: cur_vgi, cur_ind_sorted, cur_mgi

INTEGER, ALLOCATABLE                    :: counter_A(:), counter_B(:)
INTEGER                                 :: n_assoc

LOGICAL                                 :: confirmed, lower_resolution

 confirmed = .FALSE.
 lower_resolution = .FALSE.
 DO WHILE(.NOT. confirmed) ! #00 main loop
  ! getting number of virtGrid cells
  IF(lower_resolution) THEN
   IF(n_in_cell > 1) THEN
    n_in_cell = n_in_cell - 1
    lower_resolution = .FALSE.
   ELSE
    EXIT
   END IF
  ! the initial estimate
  ELSE
   n_in_cell = FLOOR((n_modelgrid/min_incell)**(1.0/3.0))
  END IF

  ! write(*,*) 'virt_gridAB_init: n_in_cell = ', n_in_cell

  IF(n_in_cell < 1) n_in_cell = 1
  
  ! we allocate arays only for modGrid cells associated to at least one propGrid cell
  ! x, y, z represent here generalized coordinates
  ! 1D: x ~ radius, y,z: no meaning
  ! 2D: x ~ radius, y ~ angle, z: no meaning
  ! 3D: x, y, z: standard meaning
  n_assoc = 0
  DO cur_mgi = 1, n_modelgrid - add_mg
   IF(model_grid(cur_mgi)%assoc_cells > 0) n_assoc = n_assoc + 1
  END DO
  ! write(*,*) 'virt_gridAB_init: n_assoc = ', n_assoc
  N_vgrid_x = n_in_cell
  N_vgrid_y = n_in_cell
  N_vgrid_z = n_in_cell
  N_vgrid_cells_A = N_vgrid_x * N_vgrid_y * N_vgrid_z
  N_vgrid_cells_B = (N_vgrid_x - 1) * (N_vgrid_y - 1) * (N_vgrid_z - 1)
  ALLOCATE(n_points_A(N_vgrid_cells_A), n_points_B(N_vgrid_cells_B + 1))
  ALLOCATE(indices_A(N_vgrid_cells_A), indices_B(N_vgrid_cells_B + 1))
  ALLOCATE(counter_A(N_vgrid_cells_A), counter_B(N_vgrid_cells_B + 1))
  ALLOCATE(vg_indexy_A(n_assoc, 2), vg_indexy_B(n_assoc, 2))
  ALLOCATE(vg_pom_A(n_assoc, 2), vg_pom_B(n_assoc, 2))
  n_points_A(:) = 0
  n_points_B(:) = 0
  n_zeros = 0
  
  IF(model_type == 1) THEN
   vg_xmin = R_star
   vg_xmax = R_inf
   w_vgrid_x = (R_inf - R_star)/N_vgrid_x
   vg_ymin = ymin
   vg_ymax = ymax
   w_vgrid_y = (R_inf - R_star)
   vg_zmin = zmin
   vg_zmax = zmax
   w_vgrid_z = (R_inf - R_star)
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
  cur_index_mgi = 0
  DO cur_point = 1, n_modelgrid - add_mg
   ! write(*,*) 'virt_gridAB_init: assoc_cells = ', model_grid(cur_point)%assoc_cells
   IF(model_grid(cur_point)%assoc_cells <= 0) CYCLE
   cur_index_mgi = cur_index_mgi + 1

   IF(model_type == 1) THEN
    cur_x = model_grid(cur_point)%rwind
    cur_y = 0.D0
    cur_z = 0.D0
    ! write(*,*) 'virt_gridAB_init: cur_x = ', cur_x
   ELSE IF(model_type == 2) THEN
    IF(inputmodel == 1) THEN
     cur_x = model_grid(cur_index_mgi)%rwind
     cur_y = model_grid(cur_point)%angle
     cur_z = 0.D0
     ! write(*,*) 'virt_gridAB_init: cur_point = ', cur_point, ' rwind = ', model_grid(cur_point)%rwind
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
    ! write(*,*) 'virt_gridAB_init: cur_n_x_A = ', cur_n_x_A
   ELSE IF(model_type == 2) THEN
    cur_n_x_A = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
    cur_n_y_A = floor((cur_y - vg_ymin)/w_vgrid_y) + 1
    cur_n_z_A = 1
    ! write(*,*) 'virGrid: cur_x/R_star = ', cur_x/R_star
   ELSE IF(model_type == 3) THEN
    cur_n_x_A = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
    cur_n_y_A = floor((cur_y - vg_ymin)/w_vgrid_y) + 1
    cur_n_z_A = floor((cur_z - vg_zmin)/w_vgrid_z) + 1
   END IF
  
   ! n_A -- numerical index of VG cell
   IF(model_type == 1) THEN
    n_A = cur_n_x_A
    IF(n_A == 0) THEN
     STOP 'virt_gridAB_init: n_A = 0'
    END IF
   ELSE IF(model_type == 2) THEN
    n_A = cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1)
   ELSE IF(model_type == 3) THEN
    n_A = cur_n_x_A + N_vgrid_x * (cur_n_y_A - 1) + N_vgrid_x * N_vgrid_y * (cur_n_z_A - 1)
   END IF
   ! write(*,*) 'virt_gridAB_init: n_x = ', cur_n_x_A, ' n_y = ', cur_n_y_A, ' n_z = ', cur_n_z_A
   ! write(*,*) 'virt_gridAB_init: n_A = ', n_A
   ! n_points -- number of points for the given cell
   n_points_A(n_A) = n_points_A(n_A) + 1
   ! vg_indexy -- list of indeces model grid --> VG index point
   ! write(*,*) 'virt_gridAB_init: n_A( ', cur_index_mgi, ') = ', n_A
   ! write(*,*) 'virt_gridAB_init: n_points_A(', cur_index_mgi, ') = ', n_points_A(n_A)
   vg_indexy_A(cur_index_mgi, 1) = n_A
   vg_indexy_A(cur_index_mgi, 2) = cur_point
   IF(n_A > N_vgrid_cells_A) THEN
    write(*,*) 'virt_gridAB_init: cur_n_x_A = ', cur_n_x_A
    write(*,*) 'virt_gridAB_init: cur_n_y_A = ', cur_n_y_A
    write(*,*) 'virt_gridAB_init: cur_n_z_A = ', cur_n_z_A
    STOP 'virt_gridAB_init n_A > N_vgrid_cells_A'
   END IF
   ! write(*,*) 'virt_gridAB_init: n_A = ', n_A, ' cur_point = ', cur_point
   !!!!!!!!
   ! repete for the B grid
   IF(cur_x > vg_xmin + w_vgrid_x/2.0 .and. cur_x < vg_xmax - w_vgrid_x/2.0 .and.&
    & cur_y > vg_ymin + w_vgrid_y/2.0 .and. cur_y < vg_ymax - w_vgrid_y/2.0 .and. &
    & cur_z > vg_zmin + w_vgrid_z/2.0 .and. cur_z < vg_zmax - w_vgrid_z/2.0 ) THEN
    IF(model_type == 1) THEN
     cur_n_x_B = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
     ! write(*,*) 'virt_gridAB_init: cur_n_x_B = ', cur_n_x_B
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
    vg_indexy_B(cur_index_mgi, 1) = n_B
    vg_indexy_B(cur_index_mgi, 2) = cur_point
   ELSE
    vg_indexy_B(cur_index_mgi, 1) = 0
    vg_indexy_B(cur_index_mgi, 2) = cur_point
    n_zeros = n_zeros + 1
   END IF
  END DO
  
  ! test of the grid
  IF(MINVAL(n_points_A, MASK=(n_points_A > 0)) < 2 * model_type .or. &
   MINVAL(n_points_B, MASK=(n_points_B > 0)) < 2 * model_type) THEN
   lower_resolution = .TRUE.
   DEALLOCATE(n_points_A, n_points_B, indices_A, indices_B)
   DEALLOCATE(counter_A, counter_B, vg_indexy_A, vg_indexy_B)
   DEALLOCATE(vg_pom_A, vg_pom_B)
  ELSE
   lower_resolution = .FALSE.
   confirmed = .TRUE.
  END IF
 END DO ! #00 end of main loop
 ! write(*,*) 'virt_gridAB_init: vg_indexy_B = ', vg_indexy_B
 
 !_______________________________________________________________
 !    #02            SORTING
 !_______________________________________________________________
 ! sort the vg_indexy according to the VG index
 ! write(*,*) 'virt_gridAB_init: started sorting A'
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

 ! write(*,*) 'virt_gridAB_init: min A: ', MINVAL(n_points_A, MASK=(n_points_A > 0)),&
 ! ' min B: ', MINVAL(n_points_B, MASK=(n_points_B > 0))

 counter_A(:) = n_points_A(:)
 counter_B(:) = n_points_B(:)
 counter_B(N_vgrid_cells_B + 1) = n_zeros
 ! write(*,*) 'virt_gridAB_init: counter_A = ', counter_A

 cur_point = 0
 DO cur_mgi = 1, n_modelgrid - add_mg
  IF(model_grid(cur_mgi)%assoc_cells <= 0) CYCLE
  cur_point = cur_point + 1
  ! A
  cur_vgi = vg_indexy_A(cur_point, 1)
  ! write(*,*) 'virt_gridAB_init: vg_indexy_A = ', cur_vgi
  IF(counter_A(cur_vgi) > 0) THEN
   cur_ind_sorted = indices_A(cur_vgi) + counter_A(cur_vgi) - 1
   vg_pom_A(cur_ind_sorted,:) = vg_indexy_A(cur_point,:)
  ELSE
   write(*,*) 'virt_gridAB_init: sorting is not OK'
   write(*,*) 'want to add a point of cur_vgi = ', cur_vgi
   write(*,*) 'with no point left'
   STOP 'virt_gridAB_init'
  END IF

  counter_A(cur_vgi) = counter_A(cur_vgi) - 1
  ! B
  cur_vgi = vg_indexy_B(cur_point, 1)
  ! write(*,*) 'virt_gridAB_init: cur_vgi = ', cur_vgi
  IF(cur_vgi > 0 ) THEN
   IF(counter_B(cur_vgi) > 0) THEN
    cur_ind_sorted = indices_B(cur_vgi) + counter_B(cur_vgi) - 1
    vg_pom_B(cur_ind_sorted,:) = vg_indexy_B(cur_point,:)
   ELSE
    write(*,*) 'virt_gridAB_init: sorting is not OK'
    write(*,*) 'want to add a point of cur_vgi = ', cur_vgi
    write(*,*) 'with no point left'
    STOP 'virt_gridAB_init'
   END IF
   counter_B(cur_vgi) = counter_B(cur_vgi) - 1
  ELSE IF(cur_vgi == 0) THEN
   cur_ind_sorted = counter_B(N_vgrid_cells_B + 1) - 1
   counter_B(N_vgrid_cells_B + 1) = counter_B(N_vgrid_cells_B + 1) - 1
  END IF
 END DO

 ! vg_indexy_A = vg_pom_A
 ! vg_indexy_B = vg_pom_B

 ! write(*,*) 'virt_gridAB_init: vg_indexy_B = ', vg_indexy_B(:,1)

 
 !_______________________________________________________________
 ! index array
 !_______________________________________________________________
 

 ! write(*,*) 'virt_gridAB_init: vg_indexy_A = ', vg_indexy_A
 is_initialized = .true.

END SUBROUTINE virt_gridAB_init

END MODULE virt_gridAB
