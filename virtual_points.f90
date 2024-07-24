! this subroutine creates virtual points for a dynamical grid calculation
SUBROUTINE virtual_points(dimIM)

 USE types
USE constants

 IMPLICIT NONE 

 ! input variables
 INTEGER                        :: dimIM
 ! 
 INTEGER                        :: I,J,NP
 ! 1D model: intervals for point distribution
 DOUBLE PRECISION               :: radius, phi, theta, angle
 DOUBLE PRECISION, DIMENSION(3) :: direction
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
 DOUBLE PRECISION, DIMENSION(3) :: pos, width

 INTEGER                        :: ind_cell_numb, Npoint
 INTEGER                        :: index_x, index_y, index_z

 TYPE(virt_point)               :: dummy
 DOUBLE PRECISION               :: A

 DOUBLE PRECISION               :: suma



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
 DO I = 1, n_modelgrid
  sumr = sumr + (model_grid(I)%rwind / R_inf) ** delta
 END DO
 ! now we will compute given numbers of points for the given spheres
 DO I = 1, n_modelgrid
  nOfPoints(I) = INT(FLOAT(Nvirtpoint) * (model_grid(I)%rwind / R_inf) ** delta / sumr)
 END DO
 DO J = 1, n_modelgrid 
  sumpart = sumpart + nOfPoints(J)
 END DO
 zbytek = Nvirtpoint - sumpart
 nOfPoints(n_modelgrid) = nOfPoints(n_modelgrid) + zbytek
 sumpart = 0
 DO J = 1, n_modelgrid 
  sumpart = sumpart + nOfPoints(J)
 END DO
 ! printing number of points for each model grid
! OPEN(UNIT=8,FILE='vp_distribution.dat')
!  DO I = 1, n_modelgrid
!   write(8,*) I, nOfPoints(I)
!  END DO
! CLOSE(8)
 ! we have zero point located
 NP = 0
 ! distribution of point on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  np_shell = nOfPoints(I)
  !IF (np_shell == 0) STOP 'number of virtual point is small'
  DO J = 1, np_shell
   NP = NP + 1
   CALL random_unitvector(direction)
   virtual_point(NP)%pos = radius * direction
   ! write(20,*) virtual_point(NP)%pos(1), virtual_point(NP)%pos(2), virtual_point(NP)%pos(3)
  END DO
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! 2D MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 delta = 1.D-1
 ! Nvirtpoint = 100 * n_modelgrid
 ALLOCATE (virtual_point(Nvirtpoint))
 write(99,*) 'number of point: ', Nvirtpoint
 write(99,*) 'computing positions of virtual point...'
 rhomax = MAXVAL(model_grid(:)%rho)
 NP = 0
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 sumr = 0.D0
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  angle = model_grid(I)%angle
  sumr = sumr + (2.0 * pi * radius * cos(angle)) ** delta
 END DO
 suma = 0
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  angle = model_grid(I)%angle
  nOfPoints(I) = FLOOR(FLOAT(Nvirtpoint) * (2.0 * pi * radius * cos(angle)) ** delta / sumr)
  suma = suma + nOfPoints(I)
  if(suma > Nvirtpoint) then
   write(*,*) 'virtual_points: I = ', I, ' z ', n_modelgrid
   STOP 'suma > Nvirtpoint'
  end if
 END DO
 ! now we will compute given numbers of points for the given spheres
 ! we have zero point located
 NP = 0
 ! distribution of point on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  np_shell = nOfPoints(I)
  !IF (np_shell == 0) STOP 'number of virtual point is small'
  DO J = 1, np_shell
   NP = NP + 1
   phi = 2.D0*pi*ran2(idum)
   theta = model_grid(I)%angle
   virtual_point(NP)%pos(1) = radius * cos(theta) * cos(phi)
   virtual_point(NP)%pos(2) = radius * cos(theta) * sin(phi)
   virtual_point(NP)%pos(3) = radius * sin(theta)
!  write(20,*) virtual_point(NP)%pos(1), virtual_point(NP)%pos(2), virtual_point(NP)%pos(3)
  END DO
 END DO
CASE(3)
 ALLOCATE(virtual_point(n_modelgrid))
 write(*,*) 'virtual_point: dim vp = ', SIZE(virtual_point)
 DO I = 1, n_modelgrid
  virtual_point(I)%pos = model_grid(I)%vec_pos
  virtual_point(I)%ind_mcell = I
 END DO
CASE DEFAULT
 STOP 'wrong choice of input model dimension...'
END SELECT
! CLOSE(20)

! every virtual point is located in the basic cell, which index can be already estimated
Npoint = SIZE(virtual_point)
DO I = 1, Npoint
 pos = virtual_point(I)%pos
 width = dyn_cell(1)%width
 index_x = FLOOR(pos(1)/width(1) + DBLE(nx_cell)/2) + 1
 index_y = FLOOR(pos(2)/width(2) + DBLE(ny_cell)/2) + 1
 index_z = FLOOR(pos(3)/width(3) + DBLE(nz_cell)/2) + 1
 ind_cell_numb = (index_x - 1) * ny_cell * nz_cell + (index_y - 1) * nz_cell + index_z
 virtual_point(I)%ind_pcell = ind_cell_numb
 dyn_cell(ind_cell_numb)%n_virt = dyn_cell(ind_cell_numb)%n_virt + 1
END DO

! sort virtual points by its number
DO J = 2, Npoint
 I = J - 1
 A = virtual_point(J)%ind_pcell
 DO WHILE (I .GE. 1)
  IF(virtual_point(I)%ind_pcell > A) THEN
   dummy = virtual_point(I + 1)
   virtual_point(I + 1) = virtual_point(I) 
   virtual_point(I) = dummy
  END IF
   I = I - 1
 END DO
END DO
! DO I = 1, Npoint
!  write(*,*) 'virtual_point: sorted vp: cell: ', virtual_point(I)%ind_pcell
! END DO




END SUBROUTINE virtual_points
