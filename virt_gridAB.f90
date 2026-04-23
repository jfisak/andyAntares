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

INTEGER, PARAMETER                      :: ind_mg = 2, ind_vg = 1

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

DOUBLE PRECISION                        :: cur_x, cur_y, cur_z

INTEGER                                 :: cur_vgi, cur_ind_sorted, cur_mgi

INTEGER, ALLOCATABLE                    :: counter_A(:), counter_B(:)
INTEGER                                 :: n_assoc

LOGICAL                                 :: confirmed, lower_resolution

DOUBLE PRECISION, DIMENSION(const_dimofspace) :: cur_pos
LOGICAL                                 :: is_vgrid_A

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
   n_in_cell = FLOOR((n_modelgrid/min_incell)**0.5)
   ! n_in_cell = 5
  END IF

  IF(n_in_cell < 1) n_in_cell = 1
  
  ! we allocate arays only for modGrid cells associated to at least one propGrid cell
  ! x, y, z represent here generalized coordinates
  ! 1D: x ~ radius, y,z: no meaning
  ! 2D: x ~ radius, y ~ angle, z: no meaning
  ! 3D: x, y, z: standard meaning
  ! n_assoc = 0
  n_assoc = n_modelgrid
  ! DO cur_mgi = 1, n_modelgrid - add_mg
  !  IF(model_grid(cur_mgi)%assoc_cells > 0) n_assoc = n_assoc + 1
  ! END DO
  N_vgrid_x = n_in_cell
  N_vgrid_y = n_in_cell
  N_vgrid_z = n_in_cell
  IF(model_type == 1) THEN
   N_vgrid_y = 1
   N_vgrid_z = 1
  ELSE IF(model_type == 2) THEN
   N_vgrid_z = 1
  END IF

  ! write(*,*) 'virt_gridAB_init: N_vgrid_x = ', N_vgrid_x, ' N_vgrid_y = ', N_vgrid_y
  ! write(*,*) 'virt_gridAB_init: allocating arrays'
  N_vgrid_cells_A = N_vgrid_x * N_vgrid_y * N_vgrid_z
  IF(model_type == 1 .AND. N_vgrid_x > 1) THEN
   N_vgrid_cells_B = N_vgrid_x - 1
  ELSE IF(model_type == 2 .AND. N_vgrid_x > 1 .AND. N_vgrid_y > 1) THEN
   N_vgrid_cells_B = (N_vgrid_x - 1) * (N_vgrid_y - 1)
  ELSE IF(model_type == 3 .AND. N_vgrid_x > 1 .AND. N_vgrid_y > 1 .AND. N_vgrid_z > 1) THEN
   N_vgrid_cells_B = (N_vgrid_x - 1) * (N_vgrid_y - 1) * (N_vgrid_z - 1)
  END IF
  ALLOCATE(n_points_A(N_vgrid_cells_A), n_points_B(N_vgrid_cells_B + 1))
  ALLOCATE(indices_A(N_vgrid_cells_A), indices_B(N_vgrid_cells_B))
  ALLOCATE(counter_A(N_vgrid_cells_A), counter_B(N_vgrid_cells_B + 1))
  ALLOCATE(vg_indexy_A(n_assoc, 2), vg_indexy_B(n_assoc, 2))
  ALLOCATE(vg_pom_A(n_assoc, 2), vg_pom_B(n_assoc, 2))
  vg_indexy_A(:,:) = 0
  vg_indexy_B(:,:) = 0
  indices_A(:) = 0
  indices_B(:) = 0
  vg_pom_A(:,:) = 0
  vg_pom_B(:,:) = 0
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
    vg_xmin = MINVAL(model_grid(:)%rwind)
    vg_xmax = R_inf
    w_vgrid_x = (R_inf - R_star)/N_vgrid_x
    vg_ymin = 0.D0
    vg_ymax = const_pi+1.D-2
    w_vgrid_y = const_pi/N_vgrid_y
    vg_zmin = -1.D0
    vg_zmax = 1.D0
    w_vgrid_z = 0.D0
   ELSE IF(inputmodel == 2) THEN
    vg_xmin = MINVAL(model_grid(:)%rxywind)
    vg_xmax = MAXVAL(model_grid(:)%rxywind)
    w_vgrid_x = (vg_xmax - vg_xmin)/N_vgrid_x
    vg_ymin = MINVAL(model_grid(:)%zwind)
    vg_ymax = MAXVAL(model_grid(:)%zwind)
    w_vgrid_y = (vg_ymax - vg_ymin)/N_vgrid_y
    vg_zmin = -1.D0
    vg_zmax = 1.D0
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
   IF(model_type == 1) THEN
    cur_x = model_grid(cur_point)%rwind
    cur_y = 0.D0
    cur_z = 0.D0
   ELSE IF(model_type == 2) THEN
    IF(inputmodel == 1) THEN
     cur_x = model_grid(cur_point)%rwind
     cur_y = model_grid(cur_point)%angle
     cur_z = 0.D0
    ELSE IF(inputmodel == 2) then
     cur_x = model_grid(cur_point)%rxywind
     cur_y = model_grid(cur_point)%zwind
     cur_z = 0.D0
    END IF
   ELSE IF(model_type == 3) THEN
    cur_x = model_grid(cur_point)%vec_pos(ind_x)
    cur_y = model_grid(cur_point)%vec_pos(ind_y)
    cur_z = model_grid(cur_point)%vec_pos(ind_z)
   END IF

   cur_pos = (/ cur_x, cur_y, cur_z /)
   is_vgrid_A = .TRUE.
   CALL get_scalar_index(cur_pos, is_vgrid_A, n_A)

   ! n_points -- number of points for the given cell
   n_points_A(n_A) = n_points_A(n_A) + 1
   ! vg_indexy -- list of indeces model grid --> VG index point
   vg_indexy_A(cur_point, ind_vg) = n_A
   IF(n_A == 0) THEN
    STOP 'virt_gridAB_init: n_A == 0'
   ELSE IF(n_A > N_vgrid_cells_A) THEN
    STOP 'virt_gridAB_init: n_A > N_vgrid_cells_A'
   END IF
   vg_indexy_A(cur_point, ind_mg) = cur_point
   IF(n_A > N_vgrid_cells_A) THEN
    STOP 'virt_gridAB_init n_A > N_vgrid_cells_A'
   END IF
   !!!!!!!!
   ! repete for the B grid
   cur_pos = (/ cur_x, cur_y, cur_z /)
   is_vgrid_A = .FALSE.
   CALL get_scalar_index(cur_pos, is_vgrid_A, n_B)
   IF(n_B > N_vgrid_cells_B) THEN
    STOP 'n_B > N_vgrid_cells_B'
   END IF
   IF(n_B > 0) THEN
    n_points_B(n_B) = n_points_B(n_B) + 1
    vg_indexy_B(cur_point, ind_vg) = n_B
    vg_indexy_B(cur_point, ind_mg) = cur_point
   ELSE
    vg_indexy_B(cur_point, ind_vg) = 0
    vg_indexy_B(cur_point, ind_mg) = cur_point
    n_zeros = n_zeros + 1
   END IF
  END DO
  n_points_B(N_vgrid_cells_B + 1) = n_zeros
  
  ! test of the grid
  IF(MINVAL(n_points_A, MASK=(n_points_A > 0)) < 2 * model_type .or. &
   MINVAL(n_points_B, MASK=(n_points_B > 0)) < 2 * model_type) THEN
   lower_resolution = .TRUE.
   ! write(*,*) 'virt_gridAB_init: lower_resolution = ', lower_resolution
   ! write(*,*) 'virt_gridAB_init: deallocating arrays'
   DEALLOCATE(n_points_A, n_points_B, indices_A, indices_B)
   DEALLOCATE(counter_A, counter_B, vg_indexy_A, vg_indexy_B)
   DEALLOCATE(vg_pom_A, vg_pom_B)
  ELSE
   lower_resolution = .FALSE.
   confirmed = .TRUE.
  END IF
 END DO ! #00 end of main loop
 ! Check if indices_A is allocated
 IF (.NOT. ALLOCATED(indices_A)) THEN
  WRITE(*,*) 'virt_gridAB_init: Error - indices_A is not allocated at the end of the subroutine.'
  STOP 'Error: indices_A not allocated'
 END IF
 !_______________________________________________________________
 !    #02            SORTING
 !_______________________________________________________________
 ! sort the vg_indexy according to the VG index
 ! calculating number of points in each vg grid cell (A)
 ! now we know where the indeces will be located, we set up an
 ! array containing the initial indeces for all vg indeces
  cur_ind_A = 0
  cur_ind_B = n_zeros
  !_______________________________________________________________
  ! create arrays with indeces pointing to an ordered list of
  ! modCell grids indeces
  DO cur_vpg_cell = 1, N_vgrid_cells_A
   indices_A(cur_vpg_cell) = cur_ind_A + 1
   cur_ind_A = cur_ind_A + n_points_A(cur_vpg_cell)
  END DO
  !_______________________________________________________________
  DO cur_vpg_cell = 1, N_vgrid_cells_B
   indices_B(cur_vpg_cell) = cur_ind_B + 1
   cur_ind_B = cur_ind_B + n_points_B(cur_vpg_cell)
  END DO
 
 
  counter_A(:) = n_points_A(:)
  counter_B(:) = n_points_B(:)
  ! write(*,*) 'virt_gridAB_init: n_zeros = ', n_zeros
  ! write(*,*) 'virt_gridAB_init: indices_B = ', indices_B
  ! write(*,*) 'virt_gridAB_init: n_points_B = ', n_points_B
  ! STOP 'virt_gridAB_init: testing'
 
  DO cur_mgi = 1, n_modelgrid
   ! A
   cur_vgi = vg_indexy_A(cur_mgi, ind_vg)
   IF(counter_A(cur_vgi) > 0) THEN
    cur_ind_sorted = indices_A(cur_vgi) + counter_A(cur_vgi) - 1
    ! IF(cur_mgi < n_modelgrid) THEN
    !  cur_ind_sorted = cur_mgi + 1
    ! ELSE
    !  cur_ind_sorted = 1
    ! END IF
    vg_pom_A(cur_ind_sorted,:) = vg_indexy_A(cur_mgi,:)
    counter_A(cur_vgi) = counter_A(cur_vgi) - 1
   ELSE
    write(*,*) 'virt_gridAB_init: sorting is not OK'
    write(*,*) 'want to add a point of cur_vgi = ', cur_vgi
    write(*,*) 'with no point left'
    STOP 'virt_gridAB_init'
   END IF
  
  
   ! B
   cur_vgi = vg_indexy_B(cur_mgi, ind_vg)
   IF(cur_vgi > 0 ) THEN
    IF(counter_B(cur_vgi) > 0) THEN
     cur_ind_sorted = indices_B(cur_vgi) + counter_B(cur_vgi) - 1
     vg_pom_B(cur_ind_sorted,:) = vg_indexy_B(cur_mgi,:)
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
 
  ! write(*,*) 'virt_gridAB_init: vg_pom_A = ', vg_pom_A(500:520,ind_vg)
  ! write(*,*) 'virt_gridAB_init: n_points_A = ', n_points_A
  ! write(*,*) 'virt_gridAB_init: minval counter_A = ', MINVAL(counter_A)
  ! STOP 'virt_gridAB_init: testing'


  vg_indexy_A(:,:) = vg_pom_A(:,:)
  vg_indexy_B(:,:) = vg_pom_B(:,:)
 ! write(*,*) 'virt_gridAB_init: vg_indexy_A = ', vg_indexy_A(1:1000,ind_vg)
 ! write(*,*) 'virt_gridAB_init: n_points_A = ', n_points_A
 DO cur_mgi = 1, n_modelgrid
  IF(vg_indexy_A(cur_mgi, ind_mg) == 0) THEN
   write(*,*) 'virt_gridAB_init: cur_mgi = ', cur_mgi
   STOP 'virt_gridAB_init: cur_mgi == 0'
  END IF
 END DO
 ! DO cur_ind_vgi = 1, n_modelgrid
 !  cur_vgi = vg_indexy_B(cur_ind_vgi, ind_vg)
 !  cur_mgi = vg_indexy_B(cur_ind_vgi, ind_mg)
 !  IF(cur_vgi > 0) THEN
 !   cur_x = model_grid(cur_mgi)%rwind
 !   cur_y = model_grid(cur_mgi)%angle
 !   cur_z = 0.D0
 !   cur_pos = (/ cur_x, cur_y, cur_z /)
 !   is_vgrid_A = .false.
 !   CALL get_scalar_index(cur_pos, is_vgrid_A, cur_scalar_mg_vgi_index)
 !   IF(cur_vgi /= cur_scalar_mg_vgi_index) THEN
 !    write(*,*) 'virt_gridAB_init: cur_scalar_mg_vgi_index = ', cur_scalar_mg_vgi_index
 !    write(*,*) 'virt_gridAB_init: cur_mgi = ', cur_mgi, ' cur_vgi = ', cur_vgi
 !    STOP 'cur_vgi /= cur_vgi_calc'
 !   END IF
 !  END IF
 ! END DO

 ! previous_index = 0
 ! DO ind_I = 1, n_modelgrid
 !  cur_vgi = vg_indexy_B(ind_I, ind_vg)
 !  cur_index_mgi = vg_indexy_B(ind_I, ind_mg)
 !  IF(cur_vgi > 0) THEN
 !   cur_x = model_grid(cur_index_mgi)%rwind
 !   cur_y = model_grid(cur_index_mgi)%angle
 !   cur_z = 0.D0
 !   cur_pos = (/ cur_x, cur_y, cur_z /)
 !   is_vgrid_A = .false.
 !   CALL get_scalar_index(cur_pos, is_vgrid_A, cur_scalar_mg_vgi_index)
 !   cur_index = cur_scalar_mg_vgi_index
 !   IF(cur_index < previous_index) THEN
 !    write(*,*) 'virt_gridAB_init: ind_I = ', ind_I
 !    write(*,*) 'virt_gridAB_init: cur_vgi = ', cur_vgi, ' ind_I = ', ind_I
 !    write(*,*) 'virt_gridAB_init: cur_index = ', cur_index, ' previous_index = ', previous_index
 !    STOP 'cur_index < previous_index'
 !   END IF
 !   previous_index = cur_index
 !  END IF
 ! END DO
 ! STOP 'virt_gridAB_init: testing'

 ! write(*,*) 'virt_gridAB_init: vg_indexy_B = ', vg_indexy_B(:,1)

 
 !_______________________________________________________________
 ! index array
 !_______________________________________________________________
 

 ! write(*,*) 'virt_gridAB_init: vg_indexy_A = ', vg_indexy_A
 is_initialized = .true.

 ! STOP 'virt_gridAB_init: testing'

END SUBROUTINE virt_gridAB_init

SUBROUTINE get_vg_index(cur_pos, out_index_A, out_index_B)

 IMPLICIT NONE

 DOUBLE PRECISION, DIMENSION(const_dimofspace)     :: cur_pos
 DOUBLE PRECISION                                  :: cur_x, cur_y, cur_z
 INTEGER, DIMENSION(const_dimofspace)              :: out_index_A, out_index_B
 INTEGER                                :: cur_n_x_A, cur_n_y_A, cur_n_z_A
 INTEGER                                :: cur_n_x_B, cur_n_y_B, cur_n_z_B

 ! write(*,*) 'get_vg_index: 
 ! write(*,*) 'get_vg_index: cur_pos = ', (cur_pos(ind_x)-vg_xmin)/w_vgrid_x, &
 !  (cur_pos(ind_y) - vg_ymin)/w_vgrid_y, cur_pos(ind_z) - vg_zmin
 cur_x = cur_pos(ind_x)
 cur_y = cur_pos(ind_y)
 cur_z = cur_pos(ind_z)

 ! write(*,*) 'get_vg_index: vg_xmin = ', vg_xmin, ' w_vgrid_x = ', w_vgrid_x
 ! the grid A
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

 ! correction for the boundaries
 IF(cur_x >= vg_xmax .and. cur_n_x_A > N_vgrid_x) THEN
  cur_n_x_A = N_vgrid_x
 END IF
 IF(cur_y <= vg_ymax .and. cur_n_y_A > N_vgrid_y) THEN
  cur_n_y_A = N_vgrid_y
 END IF
 IF(cur_z <= vg_zmax .and. cur_n_z_A > N_vgrid_z) THEN
  cur_n_z_A = N_vgrid_z
 END IF
 
 IF(cur_n_x_A > N_vgrid_x .or. cur_n_y_A > N_vgrid_y .or. cur_n_z_A > N_vgrid_z) THEN
  write(*,*) 'get_vg_index: n_x = ', cur_n_x_A, ' n_y = ', cur_n_y_A, ' n_z = ', cur_n_z_A
  STOP 'get_vg_index: index out of range'
 END IF
 out_index_A = (/ cur_n_x_A, cur_n_y_A, cur_n_z_A /)




 ! the grid B
 IF(cur_x >= vg_xmin + w_vgrid_x/2.0 .and. cur_x <= vg_xmax - w_vgrid_x/2.0 .and.&
  & cur_y >= vg_ymin + w_vgrid_y/2.0 .and. cur_y <= vg_ymax - w_vgrid_y/2.0 .and. &
  & cur_z >= vg_zmin + w_vgrid_z/2.0 .and. cur_z <= vg_zmax - w_vgrid_z/2.0 ) THEN
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
  IF(cur_x <= vg_xmax - w_vgrid_x/2.0 .and. cur_n_x_B > N_vgrid_x - 1) THEN
   cur_n_x_B = N_vgrid_x - 1
  END IF
  IF(cur_y <= vg_ymax - w_vgrid_y/2.0 .and. cur_n_y_B > N_vgrid_y - 1 .and. model_type >= 2) THEN
   cur_n_y_B = N_vgrid_y - 1
  END IF
  IF(cur_z <= vg_zmax - w_vgrid_z/2.0 .and. cur_n_z_B > N_vgrid_z - 1 .and. model_type == 3) THEN
   cur_n_z_B = N_vgrid_z - 1
  END IF
  out_index_B = (/ cur_n_x_B, cur_n_y_B, cur_n_z_B /)
 ELSE
  out_index_B = (/ 0, 0, 0 /)
 END IF
 ! write(*,*) 'get_vg_index: out_index_B = ', out_index_B
END SUBROUTINE get_vg_index

SUBROUTINE get_vg_centre_A(cur_scal_index, out_centre)

 IMPLICIT NONE
 
 INTEGER               :: cur_scal_index
 DOUBLE PRECISION, DIMENSION(const_dimofspace)      :: out_centre
 INTEGER                                            :: cur_n_x_A, cur_n_y_A, cur_n_z_A
 INTEGER                                                :: cur_zb
 

 cur_n_z_A = (cur_scal_index - 1)/(N_vgrid_x * N_vgrid_y) + 1
 cur_zb = MOD(cur_scal_index - 1,N_vgrid_x * N_vgrid_y) + 1
 cur_n_y_A = (cur_zb - 1)/N_vgrid_x + 1
 cur_zb = MOD(cur_zb - 1, N_vgrid_x) + 1
 cur_n_x_A = cur_zb
 
 out_centre = (/ w_vgrid_x * (cur_n_x_A - 0.5) + vg_xmin, w_vgrid_y * (cur_n_y_A - 0.5) + vg_ymin, &
                  & w_vgrid_z * (cur_n_z_A - 0.5) + vg_zmin /)

END SUBROUTINE get_vg_centre_A


SUBROUTINE get_vg_centre_B(cur_scal_index, out_centre)

 USE types
 IMPLICIT NONE
 
 INTEGER              :: cur_scal_index
 DOUBLE PRECISION, DIMENSION(const_dimofspace)      :: out_centre
 INTEGER                                            :: cur_n_x_B, cur_n_y_B, cur_n_z_B
 INTEGER                                                :: cur_zb

 cur_n_z_B = (cur_scal_index - 1)/((N_vgrid_x - 1) * (N_vgrid_y - 1)) + 1
 cur_zb = MOD(cur_scal_index - 1,(N_vgrid_x - 1) * (N_vgrid_y - 1)) + 1
 cur_n_y_B = (cur_zb - 1)/(N_vgrid_x - 1) + 1
 cur_zb = MOD(cur_zb - 1, (N_vgrid_x - 1)) + 1
 cur_n_x_B = cur_zb

 IF(model_type == 1) THEN
  out_centre = (/ w_vgrid_x * cur_n_x_B + vg_xmin, 0.D0, 0.D0 /)
 ELSE IF(model_type == 2) THEN
  out_centre = (/ w_vgrid_x * cur_n_x_B + vg_xmin, w_vgrid_y * cur_n_y_B + vg_ymin, 0.D0 /)
 ELSE IF(model_type == 3) THEN
  out_centre = (/ w_vgrid_x * cur_n_x_B + vg_xmin, w_vgrid_y * cur_n_y_B + vg_ymin, &
                  w_vgrid_z * cur_n_z_B + vg_zmin /)
 END IF

END SUBROUTINE get_vg_centre_B

SUBROUTINE get_scalar_index(cur_pos, is_grid_A, out_index)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_pos
LOGICAL                                                 :: is_grid_A
INTEGER                                                 :: out_index

INTEGER                                                 :: cur_n_x, cur_n_y, cur_n_z
DOUBLE PRECISION                                        :: cur_x, cur_y, cur_z

! write(*,*) 'get_scalar_index: cur_pos = ', cur_pos(ind_x)/R_star, cur_pos(ind_y), cur_pos(ind_z)
! write(*,*) 'get_scalar_index: vg_xmin = ', vg_xmin/R_star

cur_x = cur_pos(ind_x)
cur_y = cur_pos(ind_y)
cur_z = cur_pos(ind_z)

! write(*,*) 'get_scalar_index: cur_x = ', cur_x/R_star, ' cur_y = ', cur_y, ' cur_z = ', cur_z
! write(*,*) 'get_scalar_index: N_vgrid_x = ', N_vgrid_x, 'N_vgrid_y = ', N_vgrid_y, ' N_vgrid_z = ', N_vgrid_z

! write(*,*) 'get_scalar_index: *******************************************************'
! write(*,*) 'get_scalar_index: x = ', cur_x/w_vgrid_x, ' y = ', cur_y, ' z = ', cur_z
IF(is_grid_A) THEN
 IF(model_type == 1) THEN
  cur_n_x = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
  cur_n_y = 1
  cur_n_z = 1
  ! write(*,*) 'virt_gridAB_init: cur_n_x_A = ', cur_n_x_A
 ELSE IF(model_type == 2) THEN
  cur_n_x = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
  cur_n_y = floor((cur_y - vg_ymin)/w_vgrid_y) + 1
  cur_n_z = 1
 ELSE IF(model_type == 3) THEN
  cur_n_x = floor((cur_x - vg_xmin)/w_vgrid_x) + 1
  cur_n_y = floor((cur_y - vg_ymin)/w_vgrid_y) + 1
  cur_n_z = floor((cur_z - vg_zmin)/w_vgrid_z) + 1
 END IF

 ! correction for the boundaries
 IF(cur_x >= vg_xmax .and. cur_n_x > N_vgrid_x) THEN
  cur_n_x = N_vgrid_x
 END IF
 IF(cur_y <= vg_ymax .and. cur_n_y > N_vgrid_y) THEN
  cur_n_y = N_vgrid_y
 END IF
 IF(cur_z <= vg_zmax .and. cur_n_z > N_vgrid_z) THEN
  cur_n_z = N_vgrid_z
 END IF
 ! write(*,*) 'get_scalar_index: cur_x - vg_xmin = ', cur_x - vg_xmin, ' w_vgrid_x = ', w_vgrid_x
 ! write(*,*) 'get_scalar_index: cur_y - vg_ymin = ', cur_y - vg_ymin, ' w_vgrid_y = ', w_vgrid_y
 ! write(*,*) 'get_scalar_index A: cur_n_x = ', cur_n_x, ' cur_n_y = ', cur_n_y, ' cur_n_z = ', cur_n_z
 out_index = cur_n_x + N_vgrid_x * (cur_n_y - 1) + N_vgrid_x * N_vgrid_y * (cur_n_z - 1)
 IF(out_index > N_vgrid_cells_A) THEN
  write(*,*) 'get_scalar_index: cur_n_x = ', cur_n_x, ' cur_y = ', cur_y, ' cur_n_z = ', cur_n_z
  STOP 'get_scalar_index out_index > N_vgrid_cells_A'
 END IF
ELSE ! virtgrid B
 IF(cur_x >= vg_xmin + w_vgrid_x/2.0 .and. cur_x <= vg_xmax - w_vgrid_x/2.0 .and.&
  & cur_y >= vg_ymin + w_vgrid_y/2.0 .and. cur_y <= vg_ymax - w_vgrid_y/2.0 .and. &
  & cur_z >= vg_zmin + w_vgrid_z/2.0 .and. cur_z <= vg_zmax - w_vgrid_z/2.0 ) THEN
  IF(model_type == 1) THEN
   cur_n_x = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
   cur_n_y = 1
   cur_n_z = 1
  ELSE IF(model_type == 2) THEN
   cur_n_x = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
   cur_n_y = floor((cur_y - vg_ymin)/w_vgrid_y - 1.0/2.0) + 1
   cur_n_z = 1
   ! write(*,*) 'get_scalar_index: cur_n_x = ', cur_n_x, ' cur_n_y = ', cur_n_y, &
   ! ' cur_n_z = ', cur_n_z
   ! write(*,*) 'get_scalar_index I B: cur_n_x = ', cur_n_x, ' cur_n_y = ', cur_n_y, ' cur_n_z = ', cur_n_z
  ELSE IF(model_type == 3) THEN
   cur_n_x = floor((cur_x - vg_xmin)/w_vgrid_x - 1.0/2.0) + 1
   cur_n_y = floor((cur_y - vg_ymin)/w_vgrid_y - 1.0/2.0) + 1
   cur_n_z = floor((cur_z - vg_zmin)/w_vgrid_z - 1.0/2.0) + 1
  END IF

  IF(cur_x <= vg_xmax - w_vgrid_x/2.0 .and. cur_n_x > N_vgrid_x - 1) THEN
   cur_n_x = N_vgrid_x - 1
  END IF
  IF(cur_y <= vg_ymax - w_vgrid_y/2.0 .and. cur_n_y > N_vgrid_y - 1 .and. model_type >= 2) THEN
   cur_n_y = N_vgrid_y - 1
  END IF
  IF(cur_z <= vg_zmax - w_vgrid_z/2.0 .and. cur_n_z > N_vgrid_z - 1 .and. model_type == 3) THEN
   cur_n_z = N_vgrid_z - 1
  END IF
  ! write(*,*) 'get_scalar_index: cur_n_x = ', cur_n_x, ' cur_n_y = ', cur_n_y, ' cur_n_z = ', cur_n_z
  out_index = cur_n_x + (N_vgrid_x - 1) * (cur_n_y - 1) + (N_vgrid_x - 1) * (N_vgrid_y - 1) * (cur_n_z - 1)
 ELSE
  ! write(*,*) 'get_scalar_index: index is set to zero'
  out_index = 0
  cur_n_x = 0
  cur_n_y = 0
  cur_n_z = 0
 END IF

 IF(cur_n_x > N_vgrid_x - 1) THEN
  write(*,*) 'get_scalar_index: cur_n_x = ', cur_n_x, ' N_vgrid_x = ', N_vgrid_x
  STOP 'get_scalar_index: cur_n_x > N_vgrid_x'
 END IF
 IF(model_type <= 3) THEN
  IF(cur_n_y > N_vgrid_y - 1 .and. N_vgrid_y > 1) THEN
   write(*,*) 'get_scalar_index: cur_n_y = ', cur_n_y, ' N_vgrid_y = ', N_vgrid_y
   STOP 'get_scalar_index: cur_n_y > N_vgrid_y'
  END IF
 END IF
 IF(model_type == 3) THEN
  IF(cur_n_z > N_vgrid_z - 1) THEN
   write(*,*) 'get_scalar_index: cur_n_z = ', cur_n_z, ' N_vgrid_z = ', N_vgrid_z
   STOP 'get_scalar_index: cur_n_z > N_vgrid_z'
  END IF
 END IF
 ! write(*,*) 'get_scalar_index B: cur_n_x = ', cur_n_x, ' cur_n_y = ', cur_n_y, ' cur_n_z = ', cur_n_z
END IF

! write(*,*) 'get_scalar_index: is_vgrid_A = ', is_grid_A, ' out_index = ', out_index


END SUBROUTINE get_scalar_index


END MODULE virt_gridAB
