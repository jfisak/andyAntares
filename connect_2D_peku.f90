! connects propGrid and modGrid for the case of 2D petr kurfurst model
! it develops virtual grids which significantly speeds up the calculation
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE connect_2D_peku()

USE types
USE MPI
IMPLICIT NONE

INTEGER                                 :: N_vgrid_r, N_vgrid_t
INTEGER                                 :: N_vgrid_cells_A, N_vgrid_cells_B
DOUBLE PRECISION                        :: w_vgrid_r, w_vgrid_t
INTEGER                                 :: cur_point, n_rt_A, n_rt_B
INTEGER                                 :: cur_n_r_A, cur_n_t_A, cur_n_r_B, cur_n_t_B
INTEGER                                 :: cur_ind_A, cur_ind_B, cur_vpg_cell

DOUBLE PRECISION                        :: cur_r, cur_t
DOUBLE PRECISION                        :: rmax, rmin, tmax, tmin

INTEGER, DIMENSION (n_modelgrid,2)        :: vg_indexy_A, vg_indexy_B
INTEGER, PARAMETER                      :: ind_vg = 1, ind_pg = 2
INTEGER, DIMENSION(2)                   :: dummy_var_A, dummy_A, dummy_var_B, dummy_B
INTEGER, ALLOCATABLE                    :: n_points_A(:), n_points_B(:)
INTEGER, ALLOCATABLE                    :: indices_A(:), indices_B(:)

INTEGER                                 :: ind_I
INTEGER                                 :: cur_prop_cell
INTEGER                                 :: cur_vmg_A, cur_vmg_B
INTEGER                                 :: n_zeros, n_propgrid

DOUBLE PRECISION, DIMENSION(2)          :: centre_A, centre_B
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cur_pos
DOUBLE PRECISION                        :: dist_A, dist_B

INTEGER                                 :: cur_n_points, cur_start_index, cur_end_index
INTEGER, ALLOCATABLE                    :: cur_points(:)
INTEGER, DIMENSION(n_modelgrid + add_mg) :: cur_n_assoc

INTEGER                                 :: cur_mgi, cur_VG_point
DOUBLE PRECISION                        :: cur_VG_r, cur_VG_t, delta, dist, min_point

! loop variables
! variables for calculating the shortest distance between

DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

  

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek! , N_tot_zbytek

INTEGER, DIMENSION(n_propgcells)                :: cur_model_index



! #00 Set up of virtual grids
! #01 Calculate VG index of MG cells
! #02 Sorting
! #03 Connecting

! virtual grid definition
! division of virGrid A and B
!_______________________________________________________________
!   #00             SET UP OF VIRTUAL GRIDS
!_______________________________________________________________
N_vgrid_r = 10
N_vgrid_t = 10
N_vgrid_cells_A = N_vgrid_r * N_vgrid_t
N_vgrid_cells_B = (N_vgrid_r - 1) * (N_vgrid_t - 1)
ALLOCATE(n_points_A(N_vgrid_cells_A), n_points_B(N_vgrid_cells_B))
ALLOCATE(indices_A(N_vgrid_cells_A), indices_B(N_vgrid_cells_B))
n_points_A(:) = 0
n_points_B(:) = 0
n_zeros = 0

n_propgrid = SIZE(dyn_cell)

rmax = MAXVAL(model_grid(:)%rxywind) + 1e1
rmin = MINVAL(model_grid(:)%rxywind) - 1e1
tmax = MAXVAL(model_grid(:)%zwind) + 1e-2
tmin = MINVAL(model_grid(:)%zwind) - 1e-2

w_vgrid_r = abs(rmax - rmin)/N_vgrid_r
w_vgrid_t = abs(tmax - tmin)/N_vgrid_t

!_______________________________________________________________
!   #01        CALCULATE VG INDEX OF MG CELLS
!_______________________________________________________________
! calculation of the virGrid index
DO cur_point = 1, n_modelgrid
 cur_r = model_grid(cur_point)%rxywind
 cur_t = model_grid(cur_point)%zwind
 cur_n_r_A = floor((cur_r-rmin)/w_vgrid_r) + 1
 cur_n_t_A = floor((cur_t-tmin)/w_vgrid_t) + 1
 IF(cur_r == rmax) cur_n_r_A = cur_n_r_A - 1
 IF(cur_t == tmax) cur_n_t_A = cur_n_t_A - 1

 ! n_rt_A -- numerical index of VG cell
 n_rt_A = cur_n_r_A + N_vgrid_r * (cur_n_t_A - 1)
 ! n_points -- number of points for the given cell
 n_points_A(n_rt_A) = n_points_A(n_rt_A) + 1
 ! vg_indexy -- list of indeces model grid --> VG index point
 vg_indexy_A(cur_point, ind_vg) = n_rt_A
 vg_indexy_A(cur_point, ind_pg) = cur_point
 !!!!!!!!
 ! repete for the B grid
 IF(cur_r > rmin + w_vgrid_r/2.0 .and. cur_r < rmax - w_vgrid_r/2.0 .and.&
  cur_t > tmin + w_vgrid_t/2.0 .and. cur_t < tmax - w_vgrid_t/2.0) THEN
  cur_n_r_B = floor((cur_r - rmin)/w_vgrid_r - 1.0/2.0) + 1
  cur_n_t_B = floor((cur_t - tmin)/w_vgrid_t - 1.0/2.0) + 1
  n_rt_B = cur_n_r_B + (N_vgrid_r - 1) * (cur_n_t_B - 1)
  n_points_B(n_rt_B) = n_points_B(n_rt_B) + 1
  vg_indexy_B(cur_point, ind_vg) = n_rt_B
  vg_indexy_B(cur_point, ind_pg) = cur_point
 ELSE
  vg_indexy_B(cur_point, ind_vg) = 0
  vg_indexy_B(cur_point, ind_pg) = cur_point
  n_zeros = n_zeros + 1
 END IF
END DO ! going over all modGrid cels

!_______________________________________________________________
!    #02            SORTING
!_______________________________________________________________
! sort the vg_indexy according to the VG index
DO cur_point = 2, n_modelgrid
 ind_I = cur_point - 1

 dummy_var_A = vg_indexy_A(cur_point,:)

 DO WHILE(ind_I >= 1)
  IF(vg_indexy_A(ind_I,ind_vg) > dummy_var_A(ind_x)) THEN
   ! A grid
   dummy_A = vg_indexy_A(ind_I + 1,:)
   vg_indexy_A(ind_I + 1,:) = vg_indexy_A(ind_I,:)
   vg_indexy_A(ind_I,:) = dummy_A
  END IF
  ind_I = ind_I - 1
 END DO
END DO

! sort the vg_indexy according to the VG index
DO cur_point = 2, n_modelgrid
 ind_I = cur_point - 1
 dummy_var_B = vg_indexy_B(cur_point,:)
 
 DO WHILE(ind_I >= 1)
  IF(vg_indexy_B(ind_I,ind_vg) > dummy_var_B(ind_vg)) THEN
   ! B grid
   dummy_B = vg_indexy_B(ind_I + 1,:)
   vg_indexy_B(ind_I + 1,:) = vg_indexy_B(ind_I,:)
   vg_indexy_B(ind_I,:) = dummy_B
  END IF
  ind_I = ind_I - 1
 END DO
END DO
!_______________________________________________________________
! index array
cur_ind_A = 0
cur_ind_B = n_zeros

!_______________________________________________________________
! create arrays with indeces pointing to an ordered list of modCell grids indeces
DO cur_vpg_cell = 1, N_vgrid_cells_A
 indices_A(cur_vpg_cell) = cur_ind_A + 1
 cur_ind_A = cur_ind_A + n_points_A(cur_vpg_cell)
 ! write(*,*) 'connect_2D_peku: indices_A(cur_vpg_cell) = ', indices_A(cur_vpg_cell)
END DO
!_______________________________________________________________
DO cur_vpg_cell = 1, N_vgrid_cells_B
 indices_B(cur_vpg_cell) = cur_ind_B + 1
 cur_ind_B = cur_ind_B + n_points_B(cur_vpg_cell)
END DO

! write(*,*) 'connect_2D_peku: vg_indexy_A = ', vg_indexy_A(:,1)
! 
!_______________________________________________________________
!  #03              CONNECTING
!_______________________________________________________________

! going through model point one by one and calculating the closest point

#if mpi == 1
 N_single = n_propgcells/n_tasks
 N_zbytek = n_propgcells - n_tasks * N_single
 IF(my_rank <= N_zbytek - 1) THEN
  my_start = my_rank * N_single  + my_rank + 1
  my_end = (my_rank + 1) * N_single + my_rank + 1
 ELSE IF(N_zbytek == 0) THEN
  my_start = my_rank * (N_single) + 1
  my_end = my_rank * (N_single) + N_single
 ELSE IF(my_rank > N_zbytek - 1) THEN
  my_start = my_rank * N_single  + N_zbytek + 1
  my_end = (my_rank + 1) * N_single + N_zbytek
 END IF
#else
 my_start = 1
 my_end = n_propgcells
#endif
! my_start = 1
! my_end = n_modelgrid

DO cur_prop_cell = my_start, my_end
 ! choosing the grid A or B
 ! coordinates of the PG cell
 ! is the point located inside the Vgrid?
 IF(dyn_cell(cur_prop_cell)%up_cell == 0) THEN
  cur_pos = dyn_cell(cur_prop_cell)%corner
  cur_r = sqrt(cur_pos(ind_x)**2 + cur_pos(ind_y)**2)
  cur_t = cur_pos(ind_z)
  IF(cur_r >= rmin .and. cur_r <= rmax .and. &
   cur_t >= tmin .and. cur_t <= tmax) THEN

   ! A VG
   cur_n_r_A = floor((cur_r-rmin)/w_vgrid_r) + 1
   cur_n_t_A = floor((cur_t-tmin)/w_vgrid_t) + 1
   cur_vmg_A = cur_n_r_A + N_vgrid_r * (cur_n_t_A - 1)

   ! B VG
   IF(cur_r > rmin + w_vgrid_r/2.0 .and. cur_r < rmax - w_vgrid_r/2.0 .and.&
      cur_t > tmin + w_vgrid_t/2.0 .and. cur_t < tmax - w_vgrid_t/2.0) THEN
    cur_n_r_B = floor((cur_r - rmin)/w_vgrid_r - 1.0/2.0) + 1
    cur_n_t_B = floor((cur_t - tmin)/w_vgrid_t - 1.0/2.0) + 1
    cur_vmg_B = cur_n_r_B + (N_vgrid_r - 1) * (cur_n_t_B - 1)
    ! write(*,*) 'cur_n_r_B = ', cur_n_r_B, ' cur_n_t_B = ', cur_n_t_B
   ELSE
    cur_n_t_B = 0
   END IF
   
   ! distances from the centres of the VG A and B
   centre_A(ind_x) = w_vgrid_r * (cur_n_r_B - 1) + w_vgrid_r/2.0
   centre_A(ind_y) = w_vgrid_t * (cur_n_t_B - 1) + w_vgrid_t/2.0

   centre_B(ind_x) = w_vgrid_r * (cur_n_r_A - 1) + w_vgrid_r
   centre_B(ind_y) = w_vgrid_t * (cur_n_t_A - 1) + w_vgrid_t

   ! looking for the closest point
   dist_A = sqrt((cur_r-centre_A(ind_x))**2+(cur_t-centre_A(ind_y))**2)
   IF(cur_n_t_B > 0) THEN
    dist_B = sqrt((cur_r-centre_B(ind_x))**2+(cur_t-centre_B(ind_y))**2)
   ELSE IF(cur_n_t_B == 0) THEN
    dist_B = 1.D99
   END IF

   ! choosing the correct modGrid points
   IF(dist_A < dist_B) THEN
    ! A is the winner
    ! saving modGrid points to the array
    cur_n_points = n_points_A(cur_vmg_A)
    cur_start_index = indices_A(cur_vmg_A)
    cur_end_index = cur_start_index + cur_n_points - 1

    ALLOCATE(cur_points(cur_n_points))

    cur_points = vg_indexy_A(cur_start_index:cur_end_index,ind_pg)
    
   ELSE IF(dist_A >= dist_B) THEN
    ! B is the winner
    ! saving modGrid points to the array
    cur_n_points = n_points_B(cur_vmg_B)
    cur_start_index = indices_B(cur_vmg_B)
    cur_end_index = cur_start_index + cur_n_points - 1

    ALLOCATE(cur_points(cur_n_points))

    cur_points = vg_indexy_B(cur_start_index:cur_end_index,ind_pg)
   END IF

   ! finally, looking for the point with the shortest distance
   delta = 1.D99
   DO cur_VG_point = 1, cur_n_points
    cur_mgi = cur_points(cur_VG_point)
    cur_VG_r = model_grid(cur_VG_point)%rxywind
    cur_VG_t = model_grid(cur_VG_point)%zwind
    dist = sqrt((cur_VG_r - cur_r)**2 + (cur_VG_t - cur_t)**2)
    IF(dist < delta) THEN
     delta = dist
     min_point = cur_mgi
     IF(min_point < 0) THEN
      write(*,*) 'connect_2D_peku: cur_prop_cell = ', cur_prop_cell
      STOP 'min_point < 0'
     END IF
    END IF
   END DO
   ! Heureka! We have got the point!
   cur_model_index(cur_prop_cell) = INT(min_point)
   IF(INT(min_point) /= 0) THEN
    cur_n_assoc(INT(min_point)) = cur_n_assoc(INT(min_point)) + 1
   END IF

   IF(cur_r < R_star) THEN
    cur_model_index(cur_prop_cell) = photosphere_index
    cur_n_assoc(photosphere_index) = cur_n_assoc(photosphere_index) + 1
   ELSE IF(cur_r > R_inf) THEN
    cur_model_index(cur_prop_cell) = outerspace_index
    cur_n_assoc(outerspace_index) = cur_n_assoc(outerspace_index) + 1
   ELSE IF(cur_model_index(cur_prop_cell) == 0) THEN
    cur_model_index(cur_prop_cell) = vacuum_index
    cur_n_assoc(vacuum_index) = cur_n_assoc(vacuum_index) + 1
   END IF

   DEALLOCATE(cur_points)
  ! it is outside the model grid
  ELSE  
   cur_model_index(cur_prop_cell) = outerspace_index
  END IF ! if inside the modGrid area
 ELSE ! up_cell != 0
   cur_model_index(cur_prop_cell) = -1
 END IF ! up_cell == 0
 IF(cur_model_index(cur_prop_cell) < 0) THEN
  write(*,*) 'connect_2D_peku: cur_prop_cell = ', cur_prop_cell
  STOP 'model_index < 0'
 END IF
 
END DO ! loop over every propGrid cell to calculate associated modGrid cells


#if mpi == 1
IF(n_tasks > 1) THEN
 CALL MPI_ALLREDUCE(cur_model_index(:), dyn_cell(:)%model_index, n_propgcells, &
   MPI_INT, MPI_SUM, mpi_comm_world, ierr)
 CALL MPI_ALLREDUCE(cur_n_assoc(:), model_grid(:)%assoc_cells, n_modelgrid + add_mg, &
   MPI_INT, MPI_SUM, mpi_comm_world, ierr)
ELSE IF(n_tasks == 1) THEN
 dyn_cell(:)%model_index = cur_model_index(:)
 model_grid(:)%assoc_cells = cur_n_assoc(:)
END IF

#endif





END SUBROUTINE connect_2D_peku
