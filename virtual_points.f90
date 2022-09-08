! this subroutine creates virtual points for a dynamical grid calculation
SUBROUTINE virtual_points(dimIM)

 USE types

 IMPLICIT NONE 

 ! input variables
 INTEGER                        :: dimIM
 ! 
 INTEGER                        :: I,J,NP
 ! 1D model: intervals for point distribution
 DOUBLE PRECISION               :: radius, phi
 DOUBLE PRECISION, DIMENSION(3) :: direction
 INTEGER                        :: np_shell
 ! division of an interval [0, 1] into parts corresponding to a density
 DOUBLE PRECISION               :: rhomax
 ! virtual point distribution
 DOUBLE PRECISION               :: sumr
 DOUBLE PRECISION               :: delta
 ! bound of the division
 DOUBLE PRECISION, DIMENSION(n_modelgrid) :: bounds
 INTEGER, DIMENSION(n_modelgrid) :: nOfPoints
 ! a random point
 DOUBLE PRECISION               :: point
 DOUBLE PRECISION               :: ran2
 INTEGER                        :: sumpart = 0, zbytek
 DOUBLE PRECISION, DIMENSION(3) :: pos, width

 INTEGER                        :: cur_ncell, cur_nvp
 INTEGER                        :: ind_cell_numb, n_cellsa, Npoint
 INTEGER                        :: ind_x, ind_y, ind_z
 INTEGER                        :: cur_order, n_cells


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
 ALLOCATE (virtual_point(Nvirtpoint))
 write(99,*) 'number of point: ', Nvirtpoint
 write(99,*) 'computing positions of virtual point...'
 rhomax = MAXVAL(model_grid(:)%rho)
 NP = 0
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 sumr = 0.D0
 DO I = 1, n_modelgrid
  sumr = sumr + (model_grid(I)%rwind / R_inf) ** delta
 END DO
 DO I = 1, n_modelgrid
  nOfPoints(I) = INT(FLOAT(Nvirtpoint) * (model_grid(I)%rwind / R_inf) ** delta / sumr)
 END DO
 ! printing number of points for each model grid
 ! OPEN(UNIT=8,FILE='vp_distribution.dat')
  DO I = 1, n_modelgrid
   write(8,*) I, nOfPoints(I)
  END DO
 CLOSE(8)
 ! now we will compute given numbers of points for the given spheres
 DO I = 1, Nvirtpoint
  point = 1.D2 * ran2(idum)
  DO J = 1, n_modelgrid
   IF((bounds(J) > point)) THEN
    nOfPoints(J) = nOfPoints(J) + 1
    EXIT
   END IF
  END DO
 END DO
 ! we have zero point located
 NP = 0
 ! distribution of point on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = sqrt(model_grid(I)%rwind**2 - model_grid(I)%zwind**2)
  np_shell = nOfPoints(I)
  !IF (np_shell == 0) STOP 'number of virtual point is small'
  DO J = 1, np_shell
   NP = NP + 1
   phi = 2.D0*pi*ran2(idum)
   virtual_point(NP)%pos(1) = radius * cos(phi)
   virtual_point(NP)%pos(2) = radius * sin(phi)
   virtual_point(NP)%pos(3) = model_grid(I)%zwind
!    write(20,*) virtual_point(NP)%pos(1), virtual_point(NP)%pos(2), virtual_point(NP)%pos(3)
  END DO
 END DO
CASE(3)
 ALLOCATE(virtual_point(n_modelgrid))
 DO I = 1, n_modelgrid
  virtual_point(I)%pos = model_grid(I)%vec_pos
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
 ind_x = FLOOR(pos(1)/width(1) + DBLE(nx_cell)/2) + 1
 ind_y = FLOOR(pos(2)/width(2) + DBLE(ny_cell)/2) + 1
 ind_z = FLOOR(pos(3)/width(3) + DBLE(nz_cell)/2) + 1
 ind_cell_numb = (ind_x - 1) * ny_cell * nz_cell + (ind_y - 1) * nz_cell + ind_z
 virtual_point(I)%n_cell = ind_cell_numb
 dyn_cell(ind_cell_numb)%n_virt = dyn_cell(ind_cell_numb)%n_virt + 1
END DO

END SUBROUTINE virtual_points
