! Set up propagation grid cells using dynamic cells
! the sbr supports the adaptive propGrid as well
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE setup_propgrid() 


USE MPI
USE types
USE constants

IMPLICIT NONE    

! loop variables
INTEGER                                :: ind_I, ind_J, K, L, cur_xyz
INTEGER, PARAMETER                     :: Nmax = 1000000000               
! variables describing dynamic cells
INTEGER                                :: max_n_dcell, N_dyn_grid
INTEGER                                :: xp, xm, yp, ym, zp, zm
TYPE(dyn_grid_cell), ALLOCATABLE       :: pom2(:)
INTEGER                                 :: status(MPI_STATUS_SIZE)
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_corner
DOUBLE PRECISION, PARAMETER             :: minival = 1.D-2
INTEGER                                 :: cur_pgi
REAL(8)                                 :: calc_x, calc_y, calc_z
INTEGER(8)                              :: ceil_x, ceil_y, ceil_z


#if mpi==1
 IF(my_rank == 0) THEN
#endif

! Number of propagation grid cells
Ngrid = nx_cell * ny_cell * nz_cell

IF (Ngrid .GT. Nmax) THEN
   PRINT*, 'ERROR: N > Nmax', Ngrid
   STOP 
END IF

! Define length of dynamic arrays (cell and package)
IF(dyngrid == 0) THEN
 ALLOCATE (dyn_cell(INT(Ngrid)))
ELSE 
 ALLOCATE (dyn_cell(INT(Ngrid + 3.0 * Nvirtpoint)))
END IF

! Size of the basic grid cells in x,y, and z direction (now they are with
! the same size i.e. regular gred)
calc_x = 2.E0 * xmax / DBLE(nx_cell)
calc_y = 2.E0 * ymax / DBLE(ny_cell)
calc_z = 2.E0 * zmax / DBLE(nz_cell)
ceil_x = CEILING(calc_x, kind=8)
ceil_y = CEILING(calc_y, kind=8)
ceil_z = CEILING(calc_z, kind=8)
basic_cell_width(ind_x) = DBLE(ceil_x)
basic_cell_width(ind_y) = DBLE(ceil_y)
basic_cell_width(ind_z) = DBLE(ceil_z)
! write(*,*) 'setup_propgrid: xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax
! write(*,*) 'setup_propgrid: calc_x = ', calc_x, ' calc_y = ', calc_y, ' calc_z = ', calc_z
! write(*,*) 'setup_propgrid: nx_cell = ', nx_cell, ' ny_cell = ', ny_cell, ' nz_cell = ', nz_cell
! write(*,'(A,F26.6,F26.6,F26.6)') 'setup_propgrid: basic_cell_width  = ', basic_cell_width
! STOP 'setup_propgrid: testing'


L = 1
DO ind_I=1, nx_cell
 DO ind_J=1, ny_cell
  DO K=1, nz_cell
   ! Index(number) of each cell in x,y, and z direction
   !dyn_cell(L)%indexc(1) = I
   !dyn_cell(L)%indexc(2) = ind_J
   !dyn_cell(L)%indexc(3) = K
   ! Coordinates of the lower left corner of each cell
   dyn_cell(L)%corner(ind_x)  = - xmax + DBLE((ind_I - 1)) * basic_cell_width(ind_x)
   dyn_cell(L)%corner(ind_y)  = - ymax + DBLE((ind_J - 1)) * basic_cell_width(ind_y)     
   dyn_cell(L)%corner(ind_z)  = - zmax + DBLE((K - 1)) * basic_cell_width(ind_z) 
   ! cell width
   dyn_cell(L)%width(ind_x) = basic_cell_width(ind_x)
   dyn_cell(L)%width(ind_y) = basic_cell_width(ind_y)
   dyn_cell(L)%width(ind_z) = basic_cell_width(ind_z)
   ! number of down cell is equal to zero
   dyn_cell(L)%down_cell = 0
   dyn_cell(L)%up_cell = 0
   dyn_cell(L)%n_sbgr = (/ 0, 0, 0 /)
   ! calculation of neighbors
   ! x+
   IF(ind_I == nx_cell) THEN
    xp = -99
   ELSE
    xp = L + ny_cell * nz_cell
   END IF
   ! x-
   IF(ind_I == 1) THEN
    xm = -99
   ELSE
    xm = L - ny_cell * nz_cell
   END IF
   ! y+
   IF(ind_J == ny_cell) THEN
    yp = -99
   ELSE
    yp = L + nz_cell
   END IF
   ! y-
   IF(ind_J == 1) THEN
    ym = -99
   ELSE
    ym = L - nz_cell
   END IF
   ! z+
   IF(K == nz_cell) THEN
    zp = -99
   ELSE
    zp = L + 1
   END IF
   ! z-
   IF(K == 1) THEN
    zm = -99
   ELSE
    zm = L - 1
   END IF
   ! 
   dyn_cell(L)%neighbor(posx) = xp
   dyn_cell(L)%neighbor(negx) = xm
   dyn_cell(L)%neighbor(posy) = yp
   dyn_cell(L)%neighbor(negy) = ym
   dyn_cell(L)%neighbor(posz) = zp
   dyn_cell(L)%neighbor(negz) = zm
   L = L + 1
  END DO
 END DO
END DO


IF (dyngrid /= 0) then
 CALL virtual_points(model_type)
END IF

! the maximal number of cells is now equal to L
max_n_dcell = Ngrid

IF(dyngrid /= 0) THEN
 DO ind_I = 1, Ngrid
    CALL create_dynamical_grid_cells(ind_I, max_n_dcell)
 END DO

! we will resize the field dyn_cell
! because we do not want empty cells
! inside
 ALLOCATE(pom2(max_n_dcell))
  pom2(:) = dyn_cell(1:max_n_dcell)
 DEALLOCATE(dyn_cell)
 ALLOCATE(dyn_cell(max_n_dcell))
  dyn_cell(:) = pom2(:)
 DEALLOCATE(pom2)
 N_dyn_grid = max_n_dcell
END IF

! definition of a total number of propGrid cells
n_propgcells = max_n_dcell
#if mpi==1
 END IF

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 IF(my_rank == 0) THEN
  DO ind_I = 1, n_tasks - 1
   CALL MPI_SEND(n_propgcells, 1, MPI_INT, ind_I, 1, MPI_COMM_WORLD, ierr)
  END DO
  write(*,*) 'setup_propgrid: my_rank = ', my_rank, ' n_propgcells = ', n_propgcells
 ELSE
  CALL MPI_RECV(n_propgcells, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, status, ierr)
  write(*,*) 'setup_propgrid: my_rank = ', my_rank, ' n_propgcells = ', n_propgcells
  ALLOCATE(dyn_cell(n_propgcells))
 END IF


 ! correction of coordinates
 ! if the corner coordinate is equal to zero, the coordinate could be set up to a really small
 ! non-zero number, we will set those numbers to zero
 DO cur_pgi = 1, n_propgcells
  cur_corner = dyn_cell(cur_pgi)%corner
  DO ind_I = 1, const_dimofspace
   IF(abs(cur_corner(ind_I)) < minival) THEN
    write(*,*) 'setup_propgrid: setting I = ', ind_I, ' corner = ', cur_corner(ind_I), ' to zero'
    dyn_cell(cur_pgi)%corner(ind_I) = 0.D0
   END IF
  END DO
 END DO

 ! sending the physical quantities to all other processes
 ! DO cur_pgi = 1, n_propgcells
 write(*,*) 'setup_propgrid: sending the propGrid informations'
 write(*,*) 'setup_propgrid: n_propgcells = ', n_propgcells
 ! vector variables
 DO cur_xyz = 1,3
  IF(my_rank == 0) THEN
   DO ind_I = 1, n_tasks - 1
    CALL MPI_SEND(dyn_cell(:)%corner(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 1, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%width(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 2, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%vec_vel(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 3, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%n_sbgr(cur_xyz), n_propgcells, MPI_INT, ind_I, 4, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(basic_cell_width(cur_xyz), 1, MPI_DOUBLE, ind_I, 5, MPI_COMM_WORLD, ierr)
   END DO
  ELSE
   CALL MPI_RECV(dyn_cell(:)%corner(cur_xyz), n_propgcells, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
   CALL MPI_RECV(dyn_cell(:)%width(cur_xyz), n_propgcells, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, status, ierr)
   CALL MPI_RECV(dyn_cell(:)%vec_vel(cur_xyz), n_propgcells, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, status, ierr)
   CALL MPI_RECV(dyn_cell(:)%n_sbgr(cur_xyz), n_propgcells, MPI_INT, 0, 4, MPI_COMM_WORLD, status, ierr)
   CALL MPI_RECV(basic_cell_width(cur_xyz), 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, status, ierr)
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 END DO
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 ! scalar variables
 IF(my_rank == 0) THEN
  DO ind_I = 1, n_tasks - 1
   CALL MPI_SEND(dyn_cell(:)%up_cell, n_propgcells, MPI_INT, ind_I, 1, MPI_COMM_WORLD, ierr)
   CALL MPI_SEND(dyn_cell(:)%down_cell, n_propgcells, MPI_INT, ind_I, 1, MPI_COMM_WORLD, ierr)
  END DO
 ELSE
  CALL MPI_RECV(dyn_cell(:)%up_cell, n_propgcells, MPI_INT, 0, 1, MPI_COMM_WORLD, status, ierr)
  CALL MPI_RECV(dyn_cell(:)%down_cell, n_propgcells, MPI_INT, 0, 1, MPI_COMM_WORLD, status, ierr)
 END IF
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 ! another vector variable with double dimension
 DO cur_xyz = 1,6
  IF(my_rank == 0) THEN
   DO ind_I = 1, n_tasks - 1
    CALL MPI_SEND(dyn_cell(:)%neighbor(cur_xyz), n_propgcells, MPI_INT, ind_I, 1, MPI_COMM_WORLD, ierr)
   END DO
  ELSE
   CALL MPI_RECV(dyn_cell(:)%neighbor(cur_xyz), n_propgcells, MPI_INT, 0, 1, MPI_COMM_WORLD, status, ierr)
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
 END DO

#endif




END SUBROUTINE setup_propgrid
