! a sbr connecting propGrid with the modGrid
! for a 2D model with defined
!
! radius and a lateral coordinate
!
! input: NONE
! output: NONE
!
SUBROUTINE connect_2D_basic()

USE MPI
USE types
USE constants
USE virt_gridAB

IMPLICIT NONE

DOUBLE PRECISION               :: diagonal
! loop variables
INTEGER                        :: cur_propcell, best_index
! variables for calculating the shortest distance between
! propagation and model cell
DOUBLE PRECISION               :: delta, delta2
! radial and vertical distance
DOUBLE PRECISION               :: pgi_radius, pgi_theta
DOUBLE PRECISION                                :: mgi_radius, mgi_theta
DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek
INTEGER                                 :: cur_mgi
DOUBLE PRECISION                :: tot_delta = 0.D0
INTEGER                         :: n_adjonced = 0

INTEGER                         :: cur_n_assoc, ind_I, n_assoc
INTEGER, DIMENSION(const_dimofspace) :: cur_index_A, cur_index_B
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_centre_A, cur_centre_B, cur_vg_index
TYPE(modelgrid), ALLOCATABLE    :: mg_pom(:)

DOUBLE PRECISION                :: dist_A, dist_B
INTEGER                         :: n_A, n_B
INTEGER                         :: cur_n_points, cur_start_index, cur_end_index

INTEGER                                 :: cur_n_x_A, cur_n_y_A, cur_n_z_A
INTEGER                                 :: cur_n_x_B, cur_n_y_B, cur_n_z_B
INTEGER                         :: cur_mgi_assoc

LOGICAL                         :: choose_A, choose_B
INTEGER                         :: cur_mgi_index
INTEGER, ALLOCATABLE            :: cur_points(:)

DOUBLE PRECISION                :: max_dist


DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_pos, cur_pg_pos
LOGICAL                                                 :: is_vgrid_A
INTEGER                                                 :: cur_scalar_vgi_index
INTEGER                                                 :: cur_scalar_pga_vgi_index, cur_scalar_pgb_vgi_index
INTEGER                                                 :: cur_scalar_mg_vgi_index
INTEGER, DIMENSION(const_dimofspace)                    :: cur_nA, cur_nB

max_dist = 2 * sqrt(basic_cell_width(ind_x)**2+basic_cell_width(ind_y)**2+basic_cell_width(ind_z)**2)

#if mpi == 1
   N_single = n_propgcells/n_tasks
   N_zbytek = n_propgcells - n_tasks * N_single
   IF(my_rank <= N_zbytek - 1) THEN
    my_start = my_rank * (N_single + 1) + 1
    my_end = my_rank * (N_single + 1) + N_single
   ELSE IF(N_zbytek == 0) THEN
    my_start = my_rank * (N_single) + 1
    my_end = my_rank * (N_single) + N_single
   ELSE
    my_start = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + 1
    my_end = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + N_single +1
   END IF
#else
   my_start = 1
   my_end = n_propgcells
#endif

! write(*,*) 'connect_2D_basic: my_rank = ', my_rank, ' n_modelgrid = ', n_modelgrid
! write(*,*) 'connect_2D_basic: my_start = ', my_start, ' my_end = ', my_end


DO cur_propcell = my_start, my_end
 ! IF(mod(cur_propcell,10000) .EQ. 0) print*, 'associating propagation grid', cur_propcell, REAL(cur_propcell)/REAL(max_n_dcell) * 1.E2, ' % completed'
 IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
   ! Absolute radius of the propagation grid cell (midle of the cell)
   pgi_radius = SQRT((dyn_cell(cur_propcell)%corner(ind_x) + dyn_cell(cur_propcell)%width(ind_x)/2.D0)**2 + &
            (dyn_cell(cur_propcell)%corner(ind_y) + dyn_cell(cur_propcell)%width(ind_y)/2.D0)**2 + &
            (dyn_cell(cur_propcell)%corner(ind_z) + dyn_cell(cur_propcell)%width(ind_z)/2.D0)**2)
   pgi_theta = acos((dyn_cell(cur_propcell)%corner(ind_z) + dyn_cell(cur_propcell)%width(ind_z)/2.D0)/pgi_radius)
   cur_pg_pos = (/ pgi_radius, pgi_theta, 0.D0 /)

  IF(pgi_radius > R_star .AND. pgi_radius < R_inf) THEN
   is_vgrid_A = .TRUE.
   CALL get_scalar_index(cur_pg_pos, is_vgrid_A, n_A)
   is_vgrid_A = .FALSE.
   CALL get_scalar_index(cur_pg_pos, is_vgrid_A, n_B)
   CALL get_vg_centre_A(cur_pg_pos, cur_centre_A)
   CALL get_scalar_index(cur_centre_A, .TRUE., cur_scalar_vgi_index)
   ! IF(cur_scalar_vgi_index /= n_B) THEN
   !  write(*,*) 'connect_2D_basic: cur_scalar_vgi_index = ', cur_scalar_vgi_index, ' n_B = ', n_B
   !  STOP 'connect_2D_basic: cur_scalar_vgi_index /= n_B'
   ! END IF
   ! write(*,*) 'connect_2D_basic: n_A = ', n_A, ' n_B = ', n_B
   IF(n_B <= 0) THEN
    choose_A = .TRUE.
    choose_B = .FALSE.
    ! cur_n_points = n_points_A(n_A)
    ! cur_start_index = indices_A(n_A)
    ! cur_end_index = cur_start_index + cur_n_points - 1
   ! 2.) we have to choose between the vg grid A and B
   ELSE ! cur_index_B(:) == 0
    ! write(*,*) 'connect_2D_basic: cur_pg_pos = ', cur_pg_pos
    CALL get_vg_centre_B(cur_pg_pos, cur_centre_B)
    CALL get_scalar_index(cur_centre_B, .FALSE., cur_scalar_vgi_index)
    IF(cur_scalar_vgi_index /= n_B) THEN
     write(*,*) 'connect_2D_basic: cur_scalar_vgi_index = ', cur_scalar_vgi_index, ' n_B = ', n_B
     STOP 'connect_2D_basic: cur_scalar_vgi_index /= n_B'
    END IF

    CALL get_scalar_index(cur_centre_B, .FALSE., cur_scalar_vgi_index)

    dist_A = sqrt(pgi_radius**2+cur_centre_A(ind_x)**2 - &
     2.0 * pgi_radius * cur_centre_A(ind_x) * cos(pgi_theta - cur_centre_A(ind_y)))
    dist_B = sqrt(pgi_radius**2+cur_centre_B(ind_x)**2 - &
     2.0 * pgi_radius * cur_centre_B(ind_x) * cos(pgi_theta - cur_centre_B(ind_y)))

    IF(dist_A < dist_B) THEN
     choose_A = .TRUE.
     choose_B = .FALSE.
    ELSE IF(dist_B < dist_A) THEN ! dist_A < dist_B
     choose_A = .FALSE.
     choose_B = .TRUE.
    ELSE ! dist_A < dist_B
     write(*,*) 'connect_2D_basic: dist_A == dist_B'
     STOP 'connect_2D_basic'
    END IF ! dist_A < dist_B
   END IF ! cur_index_B(:) == 0

   IF(choose_A) THEN
    cur_n_points = n_points_A(n_A)
    cur_start_index = indices_A(n_A)
    cur_end_index = cur_start_index + cur_n_points - 1
    ALLOCATE(cur_points(cur_n_points))
    cur_points = vg_indexy_A(cur_start_index:cur_end_index,ind_mg)
    ! testing part
    ! cur_pos = (/mgi_radius, mgi_theta, 0.D0 /)
    ! CALL get_scalar_index(cur_pos, .TRUE., cur_scalar_vgi_index)
    ! IF(cur_scalar_vgi_index /= n_A) THEN
    !  STOP 'connect_2D_basic: cur_scalar_vgi_index /= n_A'
    ! END IF
   ELSE IF(choose_B) THEN
    cur_n_points = n_points_B(n_B)
    cur_start_index = indices_B(n_B)
    cur_end_index = cur_start_index + cur_n_points - 1
    ALLOCATE(cur_points(cur_n_points))
    cur_points = vg_indexy_B(cur_start_index:cur_end_index,ind_mg)
   END IF ! choose_A or choose_B
   
   ! write(*,*) 'connect_2D_basic: cur_n_x_A = ', cur_n_x_A, ' cur_n_y_A = ', cur_n_y_A, ' cur_n_z_A = ', cur_n_z_A
   ! write(*,*) 'connect_2D_basic: -bound, point, +bound = ', vg_xmin + (cur_n_x_A - 1) * w_vgrid_x, &
   !  pgi_radius, vg_xmin + cur_n_x_A * w_vgrid_x
   ! write(*,*) 'connect_2D_basic: -bound, point, +bound = ', vg_ymin + (cur_n_y_A - 1) * w_vgrid_y, &
   !  pgi_theta, vg_ymin + cur_n_y_A * w_vgrid_y
   delta = large_number
   ! write(*,*) 'connect_2D_basic: cur_points = ', cur_points
   DO ind_I = 1, cur_n_points
    cur_mgi = cur_points(ind_I)
    mgi_radius = model_grid(cur_mgi)%rwind
    mgi_theta = model_grid(cur_mgi)%angle

    delta2 = sqrt(pgi_radius**2 + mgi_radius**2 - &
     2.0 * pgi_radius * mgi_radius * cos(pgi_theta - mgi_theta))
    IF(delta2 < delta) THEN
     delta = delta2
     cur_mgi_index = cur_mgi
     write(75,*) ind_I, abs(mgi_radius-pgi_radius)/w_vgrid_x, abs(mgi_theta-pgi_theta)/w_vgrid_y
     ! write(*,*) 'connect_2D_basic: ', ind_I, abs(mgi_radius-pgi_radius)/w_vgrid_x, &
     !  abs(mgi_theta-pgi_theta)/w_vgrid_y, delta2/R_star

    END IF
   END DO
   ! STOP 'connect_2D_basic: testing'
   ! write(*,*) 'connect_2D_basic: another point'
   ! IF(delta < max_dist) THEN
   !  mgi_radius = model_grid(cur_mgi_index)%rwind
   !  mgi_theta = model_grid(cur_mgi_index)%angle
    dyn_cell(cur_propcell)%model_index = cur_mgi_index
    model_grid(cur_mgi_index)%assoc_cells = model_grid(cur_mgi_index)%assoc_cells + 1
   ! write(*,*) 'connect_2D_basic: cur_mgi_index = ', cur_mgi_index
   ! ELSE
   !  write(*,*) 'connect_2D_basic: cur_propcell = ', cur_propcell, ' was set to be vacuum'
   !  dyn_cell(cur_propcell)%model_index = vacuum_index
   !  model_grid(vacuum_index)%assoc_cells = model_grid(vacuum_index)%assoc_cells + 1
   ! END IF
   DEALLOCATE(cur_points)

  ELSE IF(pgi_radius < R_star) THEN
   ! Cells with radius smaller than the stellar radius or larger
   ! than the winds outer radius have no associated model grid cell
   ! Make them point to the dummy model grid cell
   dyn_cell(cur_propcell)%model_index = photosphere_index     
   model_grid(photosphere_index)%assoc_cells = model_grid(photosphere_index)%assoc_cells + 1
  ELSE IF(pgi_radius > R_inf) THEN
   dyn_cell(cur_propcell)%model_index = outerspace_index     
   model_grid(outerspace_index)%assoc_cells = model_grid(outerspace_index)%assoc_cells + 1
  END IF
 END IF ! dyn_cell(cur_propcell)% up_cell == 0
 ! write(*,*) 'connect_2D_basic: cur_propcell = ', cur_propcell
 ! write(*,*) 'connect_2D_basic: modGrid index = ', dyn_cell(cur_propcell)%model_index
END DO ! loop over propGrid cells
write(99,*) 'number of propagation cells in vacuum: ', model_grid(vacuum_index)%assoc_cells
! write(*,*) 'number of propagation cells in vacuum: ', model_grid(vacuum_index)%assoc_cells, ' n_modelgrid = ', n_modelgrid, &
! ' fraction = ', model_grid(vacuum_index)%assoc_cells/REAL(n_modelgrid)
DO ind_I = 1, n_propgcells
 IF(dyn_cell(ind_I)%up_cell == 0) THEN
  cur_mgi_index = dyn_cell(ind_I)%model_index
  mgi_radius = model_grid(cur_mgi_index)%rwind
  mgi_theta = model_grid(cur_mgi_index)%angle
  pgi_radius = SQRT((dyn_cell(ind_I)%corner(ind_x) + dyn_cell(ind_I)%width(ind_x)/2.D0)**2 + &
      (dyn_cell(ind_I)%corner(ind_y) + dyn_cell(ind_I)%width(ind_y)/2.D0)**2 + &
      (dyn_cell(ind_I)%corner(ind_z) + dyn_cell(ind_I)%width(ind_z)/2.D0)**2)
  pgi_theta = acos((dyn_cell(ind_I)%corner(ind_z) + dyn_cell(ind_I)%width(ind_z)/2.D0)/pgi_radius)
  delta2 = sqrt(pgi_radius**2 + mgi_radius**2 - &
   2.0 * pgi_radius * mgi_radius * cos(pgi_theta - mgi_theta))
  write(74,*) ind_I, abs(mgi_radius-pgi_radius)/w_vgrid_x, abs(mgi_theta-pgi_theta)/w_vgrid_y, delta2
 END IF
END DO

cur_n_assoc = 0
DO ind_I = 1, n_modelgrid
 IF(model_grid(ind_I)%assoc_cells > 0) THEN
  cur_n_assoc = cur_n_assoc + model_grid(ind_I)%assoc_cells
 END IF
END DO
write(*,*) 'connect_2D_basic: cur_mgi_assoc = ', cur_mgi_assoc

IF(n_tasks > 1) THEN
 CALL MPI_ALLREDUCE(dyn_cell(:)%model_index, dyn_cell(:)%model_index, n_propgcells, &
   MPI_INT, MPI_SUM, mpi_comm_world, ierr)
 CALL MPI_ALLREDUCE(model_grid(:)%assoc_cells, model_grid(:)%assoc_cells, n_modelgrid + add_mg, &
   MPI_INT, MPI_SUM, mpi_comm_world, ierr)
END IF

DO cur_propcell = 1, n_propgcells
 IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
  IF(dyn_cell(cur_propcell)%model_index == 0) THEN
   write(*,*) 'connect_2D_basic: error, cur_propcell = ', cur_propcell
   STOP 'connect_2D_basic: model_index = 0'
  END IF
 END IF
END DO

! IF(clean_modelGrid) THEN
!  n_assoc = 0
!  ! reduce number of modGrid cells
!  DO cur_mgi = 1, n_modelgrid
!   cur_n_assoc = model_grid(cur_mgi)%assoc_cells
!   IF(cur_n_assoc > 0) THEN
!    n_assoc = n_assoc + 1
!   END IF
!  END DO
!  ALLOCATE(mg_pom(n_assoc + add_mg))
!  ind_I = 1
!  DO cur_mgi = 1, n_modelgrid + add_mg
!   cur_n_assoc = model_grid(cur_mgi)%assoc_cells
!   IF(cur_n_assoc > 0 .or. cur_mgi > n_modelgrid) THEN
!    mg_pom(ind_I) = model_grid(cur_mgi)
!    ind_I = ind_I + 1
!   END IF
!  END DO
!  DEALLOCATE(model_grid)
!  ALLOCATE(model_grid(n_assoc + add_mg))
!  model_grid = mg_pom
!  n_modelgrid = n_assoc
! END IF

! STOP 'connect_2D_basic: testing'

END SUBROUTINE connect_2D_basic
