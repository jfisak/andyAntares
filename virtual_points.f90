! this subroutine creates virtual points for a dynamical grid calculation
! 
! INPUT: dimIM(INT) -- dimension of the model grid
! OUTPUT: NONE
!
SUBROUTINE virtual_points(dimIM)

 USE types
USE constants

 IMPLICIT NONE 

 ! input variables
 INTEGER                        :: dimIM
 ! 
 INTEGER                        :: ind_I,ind_J,num_points, cur_pgi
 ! 1D model: intervals for point distribution
 DOUBLE PRECISION               :: radius, phi, theta, angle
 DOUBLE PRECISION, DIMENSION(const_dimofspace) :: direction
 INTEGER                        :: np_shell
 ! division of an interval [0, 1] into parts corresponding to a density
 DOUBLE PRECISION               :: rhomax
 ! virtual point distribution
 DOUBLE PRECISION               :: sumr
 DOUBLE PRECISION               :: delta
 INTEGER, DIMENSION(n_modelgrid) :: nOfPoints
 ! a random point
 DOUBLE PRECISION               :: ran2
 INTEGER                        :: sumpart = 0, zbytek
 DOUBLE PRECISION, DIMENSION(const_dimofspace) :: pos, width

 INTEGER                        :: ind_cell_numb, Npoint
 INTEGER                        :: index_x, index_y, index_z

 TYPE(virt_point)               :: dummy
 DOUBLE PRECISION               :: vari_A

 DOUBLE PRECISION               :: suma
 INTEGER                        :: n_virt_point
 INTEGER, ALLOCATABLE           :: counter(:), indices(:), vp_pom(:)
 INTEGER                        :: cur_ind_I, cur_ind_sorted, cur_vpi
 INTEGER                        :: N_basic_cell

 
 TYPE(virt_point), ALLOCATABLE :: pom_virtpoint(:)


! OPEN(20,FILE="virtual_point.dat")
SELECT CASE (dimIM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! 1D MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the case of 1D model
! we consider radial symmetric model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
 delta = 2.D0
 ! angle between point
 ! angleParam = 20
 ! number of virtual point
 ! Ntheta =  angleParam
 ! Nphi = 2* angleParam
 !Npart = n_modelgrid * (Ntheta - 2) * Nphi + 2 * n_modelgrid
 write(99,*) 'number of virtual point: ', Nvirtpoint
 write(99,*) 'computing positions of virtual point...'
 ALLOCATE (virtual_point(Nvirtpoint))
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 sumr = 0.D0
 DO ind_I = 1, n_modelgrid
  sumr = sumr + (model_grid(ind_I)%rwind / R_inf) ** delta
 END DO
 ! now we will compute given numbers of points for the given spheres
 DO ind_I = 1, n_modelgrid
  nOfPoints(ind_I) = INT(FLOAT(Nvirtpoint) * (model_grid(ind_I)%rwind / R_inf) ** delta / sumr)
 END DO
 DO ind_J = 1, n_modelgrid 
  sumpart = sumpart + nOfPoints(ind_J)
 END DO
 zbytek = Nvirtpoint - sumpart
 nOfPoints(n_modelgrid) = nOfPoints(n_modelgrid) + zbytek
 sumpart = 0
 DO ind_J = 1, n_modelgrid 
  sumpart = sumpart + nOfPoints(ind_J)
 END DO
 ! printing number of points for each model grid
! OPEN(UNIT=8,FILE='vp_distribution.dat')
!  DO I = 1, n_modelgrid
!   write(8,*) I, nOfPoints(I)
!  END DO
! CLOSE(8)
 ! we have zero point located
 num_points = 0
 ! distribution of point on the shell of the radius R
 DO ind_I = 1, n_modelgrid
  radius = model_grid(ind_I)%rwind
  np_shell = nOfPoints(ind_I)
  !IF (np_shell == 0) STOP 'number of virtual point is small'
  DO ind_J = 1, np_shell
   num_points = num_points + 1
   CALL random_unitvector(direction)
   virtual_point(num_points)%pos = radius * direction
   ! write(20,*) virtual_point(num_points)%pos(1), virtual_point(num_points)%pos(2), virtual_point(num_points)%pos(3)
  END DO
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! 2D MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 delta = 1.D-1
 ! Nvirtpoint = 100 * n_modelgrid
 write(99,*) 'number of point: ', Nvirtpoint
 write(99,*) 'computing positions of virtual point...'
 rhomax = MAXVAL(model_grid(:)%rho)
 num_points = 0
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 sumr = 0.D0
 ! the total "length" of 2D grid
 DO ind_I = 1, n_modelgrid
  radius = model_grid(ind_I)%rwind
  angle = model_grid(ind_I)%angle
  sumr = sumr + 2.0 * const_pi * radius * sin(angle)
 END DO

 n_virt_point = INT(sumr / (1.D1 * basic_cell_width(ind_x)))
 ! write(*,*) 'virtual_points: n_virt_point = ', n_virt_point

 suma = 0
 DO ind_I = 1, n_modelgrid
  radius = model_grid(ind_I)%rwind
  angle = model_grid(ind_I)%angle
  nOfPoints(ind_I) = CEILING(FLOAT(Nvirtpoint) * (2.0 * const_pi * radius * sin(angle)) ** delta / sumr)
  suma = suma + nOfPoints(ind_I)
 END DO
 n_virt_point = suma
 ALLOCATE(virtual_point(n_virt_point))
 ! now we will compute given numbers of points for the given spheres
 ! we have zero point located
 num_points = 0
 ! distribution of point on the shell of the radius R
 DO ind_I = 1, n_modelgrid
  radius = model_grid(ind_I)%rwind
  np_shell = nOfPoints(ind_I)
  !IF (np_shell == 0) STOP 'number of virtual point is small'
  DO ind_J = 1, np_shell
   num_points = num_points + 1
   phi = 2.D0*const_pi*ran2(idum)
   theta = model_grid(ind_I)%angle
   virtual_point(num_points)%pos(ind_x) = radius * sin(theta) * cos(phi)
   virtual_point(num_points)%pos(ind_y) = radius * sin(theta) * sin(phi)
   virtual_point(num_points)%pos(ind_z) = radius * cos(theta)
!  write(20,*) virtual_point(num_points)%pos(1), virtual_point(num_points)%pos(2), virtual_point(num_points)%pos(3)
  END DO
 END DO
CASE(3)
 ALLOCATE(virtual_point(n_modelgrid))
 write(*,*) 'virtual_point: dim vp = ', SIZE(virtual_point)
 DO ind_I = 1, n_modelgrid
  virtual_point(ind_I)%pos = model_grid(ind_I)%vec_pos
  virtual_point(ind_I)%ind_mcell = ind_I
 END DO
CASE DEFAULT
 STOP 'wrong choice of input model dimension...'
END SELECT
! CLOSE(20)

! every virtual point is located in the basic cell, which index can be already estimated
Npoint = SIZE(virtual_point)
DO ind_I = 1, Npoint
 pos = virtual_point(ind_I)%pos
 width = basic_cell_width(:)
 index_x = FLOOR(pos(ind_x)/width(ind_x) + DBLE(nx_cell)/2) + 1
 index_y = FLOOR(pos(ind_y)/width(ind_y) + DBLE(ny_cell)/2) + 1
 index_z = FLOOR(pos(ind_z)/width(ind_z) + DBLE(nz_cell)/2) + 1
 ind_cell_numb = (index_x - 1) * ny_cell * nz_cell + (index_y - 1) * nz_cell + index_z
 virtual_point(ind_I)%ind_pcell = ind_cell_numb
 dyn_cell(ind_cell_numb)%n_virt = dyn_cell(ind_cell_numb)%n_virt + 1
END DO

N_basic_cell = nx_cell * ny_cell * nz_cell

ALLOCATE(counter(N_basic_cell), indices(N_basic_cell))
ALLOCATE(vp_pom(Npoint))

counter(:) = dyn_cell(1:N_basic_cell)%n_virt

! indices of list of virtual points for each basic propGrid
cur_ind_I = 0
DO ind_I = 1, N_basic_cell
 indices(ind_I) = cur_ind_I + 1
 cur_ind_I = cur_ind_I + dyn_cell(ind_I)%n_virt
END DO


DO cur_vpi = 1, Npoint
 cur_pgi = virtual_point(cur_vpi)%ind_pcell
 IF(counter(cur_pgi) > 0) THEN
  cur_ind_sorted = indices(cur_pgi) + counter(cur_pgi) - 1
  vp_pom(cur_ind_sorted) = cur_vpi
  counter(cur_pgi) = counter(cur_pgi) - 1
 ELSE
  write(*,*) 'virt_gridAB_init: sorting is not OK'
  write(*,*) 'want to add a point of cur_vpi = ', cur_vpi
  write(*,*) 'with no point left'
  STOP 'virt_gridAB_init'
 END IF
END DO

ALLOCATE(pom_virtpoint(Npoint))
pom_virtpoint(:) = virtual_point(:)

DO ind_I = 1, Npoint
 cur_ind_sorted = vp_pom(ind_I)
 pom_virtpoint(ind_I) = virtual_point(cur_ind_sorted)
END DO

virtual_point(:) = pom_virtpoint(:)



END SUBROUTINE virtual_points
