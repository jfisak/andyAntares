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
INTEGER                                :: ind_I, ind_J, ind_K, ind_L, cur_xyz
INTEGER, PARAMETER                     :: Nmax = 1000000000               
! variables describing dynamic cells
INTEGER                                :: max_n_dcell, N_dyn_grid
INTEGER                                :: xp, xm, yp, ym, zp, zm
TYPE(dyn_grid_cell), ALLOCATABLE       :: pom2(:)
INTEGER                                 :: status(MPI_STATUS_SIZE)
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_corner
DOUBLE PRECISION, PARAMETER             :: minival = 1.D-2
DOUBLE PRECISION, PARAMETER             :: mininum = 1.D0
INTEGER                                 :: cur_pgi
REAL(8)                                 :: calc_x, calc_y, calc_z
INTEGER(8)                              :: ceil_x, ceil_y, ceil_z
LOGICAL                                 :: nx_even = .false., ny_even = .false., nz_even = .false.


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
IF(xmax < 1.D18 .and. ymax < 1.D18 .and. zmax < 1.D18) THEN
 calc_x = 2.E0 * xmax / DBLE(nx_cell)
 calc_y = 2.E0 * ymax / DBLE(ny_cell)
 calc_z = 2.E0 * zmax / DBLE(nz_cell)
 ceil_x = CEILING(calc_x, kind=8)
 ceil_y = CEILING(calc_y, kind=8)
 ceil_z = CEILING(calc_z, kind=8)
 basic_cell_width(ind_x) = DBLE(ceil_x)
 basic_cell_width(ind_y) = DBLE(ceil_y)
 basic_cell_width(ind_z) = DBLE(ceil_z)
ELSE
 basic_cell_width(ind_x) = 2.E0 * xmax / DBLE(nx_cell)
 basic_cell_width(ind_y) = 2.E0 * ymax / DBLE(ny_cell)
 basic_cell_width(ind_z) = 2.E0 * zmax / DBLE(nz_cell)
END IF
! write(*,*) 'setup_propgrid: xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax
! write(*,*) 'setup_propgrid: calc_x = ', calc_x, ' calc_y = ', calc_y, ' calc_z = ', calc_z
! write(*,*) 'setup_propgrid: nx_cell = ', nx_cell, ' ny_cell = ', ny_cell, ' nz_cell = ', nz_cell
! write(*,'(A,F26.6,F26.6,F26.6)') 'setup_propgrid: basic_cell_width  = ', basic_cell_width
! write(*,*) 'setup_propgrid: basic_cell_width  = ', basic_cell_width
! STOP 'setup_propgrid: testing'

IF(MOD(nx_cell, 2) == 0) THEN
 nx_even = .true.
END IF
IF(MOD(ny_cell, 2) == 0) THEN
 ny_even = .true.
END IF
IF(MOD(nz_cell, 2) == 0) THEN
 nz_even = .true.
END IF


ind_L = 1
DO ind_I=1, nx_cell
 DO ind_J=1, ny_cell
  DO ind_K=1, nz_cell
   ! Index(number) of each cell in x,y, and z direction
   !dyn_cell(L)%indexc(1) = I
   !dyn_cell(L)%indexc(2) = ind_J
   !dyn_cell(L)%indexc(3) = K
   ! Coordinates of the lower left corner of each cell
   dyn_cell(ind_L)%corner(ind_x) = - xmax + DBLE((ind_I - 1)) * basic_cell_width(ind_x)
   dyn_cell(ind_L)%corner(ind_y) = - ymax + DBLE((ind_J - 1)) * basic_cell_width(ind_y)     
   dyn_cell(ind_L)%corner(ind_z) = - zmax + DBLE((ind_K - 1)) * basic_cell_width(ind_z) 

   ! up corner
   dyn_cell(ind_L)%upcorner(ind_x) = -xmax + DBLE(ind_I) * basic_cell_width(ind_x)
   dyn_cell(ind_L)%upcorner(ind_y) = -ymax + DBLE(ind_J) * basic_cell_width(ind_y)     
   dyn_cell(ind_L)%upcorner(ind_z) = -zmax + DBLE(ind_K) * basic_cell_width(ind_z) 
   IF(nx_even) THEN
    IF(ind_I == nx_cell/2) THEN
     dyn_cell(ind_L)%upcorner(ind_x) = 0.D0
    END IF
    IF(ind_I == nx_cell/2 + 1) THEN
     dyn_cell(ind_L)%corner(ind_x) = 0.D0
    END IF
   END IF
   IF(ny_even) THEN
    IF(ind_J == ny_cell/2) THEN
     dyn_cell(ind_L)%upcorner(ind_y) = 0.D0
    END IF
    IF(ind_J == ny_cell/2 + 1) THEN
     dyn_cell(ind_L)%corner(ind_y) = 0.D0
    END IF
   END IF
   IF(nz_even) THEN
    IF(ind_K == nz_cell/2) THEN
     dyn_cell(ind_L)%upcorner(ind_z) = 0.D0
    END IF
    IF(ind_K == nz_cell/2 + 1) THEN
     dyn_cell(ind_L)%corner(ind_z) = 0.D0
    END IF
   END IF
   IF(dyn_cell(ind_L)%corner(ind_x) == dyn_cell(ind_L)%upcorner(ind_x) .or. &
    dyn_cell(ind_L)%corner(ind_y) == dyn_cell(ind_L)%upcorner(ind_y) .or. &
    dyn_cell(ind_L)%corner(ind_y) == dyn_cell(ind_L)%upcorner(ind_y)) THEN
    write(*,*) 'setup_propgrid: corner == upcorner'
    STOP 'setup_propgrid'
   END IF

   IF(ind_I == nx_cell) THEN
    dyn_cell(ind_L)%upcorner(ind_x) = xmax
   END IF
   IF(ind_J == ny_cell) THEN
    dyn_cell(ind_L)%upcorner(ind_y) = ymax
   END IF
   IF(ind_K == nz_cell) THEN
    dyn_cell(ind_L)%upcorner(ind_z) = zmax
   END IF

   ! write(*,*) 'setup_propgrid: corner(', ind_L, ') = ', dyn_cell(L)%corner
   ! cell width
   dyn_cell(ind_L)%width(ind_x) = basic_cell_width(ind_x)
   dyn_cell(ind_L)%width(ind_y) = basic_cell_width(ind_y)
   dyn_cell(ind_L)%width(ind_z) = basic_cell_width(ind_z)
   IF(ind_I == nx_cell .and. nx_cell > 1) THEN
    dyn_cell(ind_L)%width(ind_x) = xmax - dyn_cell(ind_L)%corner(ind_x)
   END IF
   IF(ind_J == ny_cell .and. ny_cell > 1) THEN
    dyn_cell(ind_L)%width(ind_y) = ymax - dyn_cell(ind_L)%corner(ind_y)
   END IF
   IF(ind_K == nz_cell .and. nz_cell > 1) THEN
    dyn_cell(ind_L)%width(ind_z) = zmax - dyn_cell(ind_L)%corner(ind_z)
   END IF
   ! number of down cell is equal to zero
   dyn_cell(ind_L)%down_cell = 0
   dyn_cell(ind_L)%up_cell = 0
   dyn_cell(ind_L)%n_sbgr = (/ 0, 0, 0 /)
   ! calculation of neighbors
   ! x+
   IF(ind_I == nx_cell) THEN
    xp = -99
   ELSE
    xp = ind_L + ny_cell * nz_cell
   END IF
   ! x-
   IF(ind_I == 1) THEN
    xm = -99
   ELSE
    xm = ind_L - ny_cell * nz_cell
   END IF
   ! y+
   IF(ind_J == ny_cell) THEN
    yp = -99
   ELSE
    yp = ind_L + nz_cell
   END IF
   ! y-
   IF(ind_J == 1) THEN
    ym = -99
   ELSE
    ym = ind_L - nz_cell
   END IF
   ! z+
   IF(ind_K == nz_cell) THEN
    zp = -99
   ELSE
    zp = ind_L + 1
   END IF
   ! z-
   IF(ind_K == 1) THEN
    zm = -99
   ELSE
    zm = ind_L - 1
   END IF
   ! 
   dyn_cell(ind_L)%neighbor(posx) = xp
   dyn_cell(ind_L)%neighbor(negx) = xm
   dyn_cell(ind_L)%neighbor(posy) = yp
   dyn_cell(ind_L)%neighbor(negy) = ym
   dyn_cell(ind_L)%neighbor(posz) = zp
   dyn_cell(ind_L)%neighbor(negz) = zm
   ind_L = ind_L + 1
  END DO
 END DO
END DO


IF (dyngrid /= 0) then
 CALL virtual_points(model_type)
END IF

! the maximal number of cells is now equal to ind_L
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
  ! write(*,*) 'setup_propgrid: my_rank = ', my_rank, ' n_propgcells = ', n_propgcells
 ELSE
  CALL MPI_RECV(n_propgcells, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, status, ierr)
  ! write(*,*) 'setup_propgrid: my_rank = ', my_rank, ' n_propgcells = ', n_propgcells
  ALLOCATE(dyn_cell(n_propgcells))
 END IF


 ! correction of coordinates
 ! if the corner coordinate is equal to zero, the coordinate could be set up to a really small
 ! non-zero number, we will set those numbers to zero
 DO cur_pgi = 1, n_propgcells
  cur_corner = dyn_cell(cur_pgi)%corner
  DO ind_I = 1, const_dimofspace
   IF(abs(cur_corner(ind_I)) < minival) THEN
    dyn_cell(cur_pgi)%corner(ind_I) = 0.D0
   END IF
  END DO
 END DO

 ! sending the physical quantities to all other processes
 ! DO cur_pgi = 1, n_propgcells
 ! write(*,*) 'setup_propgrid: sending the propGrid informations'
 ! write(*,*) 'setup_propgrid: n_propgcells = ', n_propgcells
 ! vector variables
 DO cur_xyz = 1,3
  IF(my_rank == 0) THEN
   DO ind_I = 1, n_tasks - 1
    CALL MPI_SEND(dyn_cell(:)%corner(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 1, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%upcorner(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 1, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%width(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 2, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%vec_vel(cur_xyz), n_propgcells, MPI_DOUBLE, ind_I, 3, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(dyn_cell(:)%n_sbgr(cur_xyz), n_propgcells, MPI_INT, ind_I, 4, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(basic_cell_width(cur_xyz), 1, MPI_DOUBLE, ind_I, 5, MPI_COMM_WORLD, ierr)
   END DO
  ELSE
   CALL MPI_RECV(dyn_cell(:)%corner(cur_xyz), n_propgcells, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
   CALL MPI_RECV(dyn_cell(:)%upcorner(cur_xyz), n_propgcells, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, status, ierr)
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

! STOP 'setup_propgrid: testing'


END SUBROUTINE setup_propgrid
