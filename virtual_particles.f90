! this subroutine creates virtual particles for a dynamical grid calculation
SUBROUTINE virtual_particles(dimIM)

 USE types

 IMPLICIT NONE 

 ! input variables
 INTEGER                        :: dimIM
 ! 
 INTEGER                        :: I,J,K,NP
 ! 1D model: intervals for particles distribution
 INTEGER                        :: Ntheta, Nphi
 DOUBLE PRECISION               :: radius, phi, theta
 DOUBLE PRECISION               :: angle
 DOUBLE PRECISION, DIMENSION(3) :: direction
 INTEGER                        :: angleParam, np_shell
 ! division of an interval [0, 1] into parts corresponding to a density
 DOUBLE PRECISION               :: rhotot, rhomax
 ! virtual particles distribution
 DOUBLE PRECISION               :: sumr
 DOUBLE PRECISION, PARAMETER             :: delta = 1.0
 ! bound of the division
 DOUBLE PRECISION               :: bound, actbound
 DOUBLE PRECISION, DIMENSION(n_modelgrid) :: bounds
 INTEGER, DIMENSION(n_modelgrid) :: nOfPoints
 ! a random point
 DOUBLE PRECISION               :: point
 DOUBLE PRECISION               :: ran2, rand_rho, rand_r
 INTEGER                        :: loc_rad


 OPEN(20,FILE="virtual_particles.dat")
SELECT CASE (dimIM)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! 1D MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the case of 1D model
! we consider radial symmetric model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
 ! angle between particles
 ! angleParam = 20
 ! number of virtual particles
 ! Ntheta =  angleParam
 ! Nphi = 2* angleParam
 !Npart = n_modelgrid * (Ntheta - 2) * Nphi + 2 * n_modelgrid
 print*, 'number of particles: ', Nvirtpart
 print*, 'computing positions of virtual particles...'
 ALLOCATE (virtual_particle(Nvirtpart))
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 sumr = 0.D0
 DO I = 1, n_modelgrid
  sumr = sumr + (model_grid(I)%rwind / R_inf) ** delta
 END DO
! ! now we divide an interval [0, 1] into parts which lenght
! ! corresponds to the density magnitude
! bound = 0.D0
! DO I = 1, n_modelgrid
!  actbound = bound - log10(model_grid(I)%rho) / rhotot
!  bounds(I) = actbound
!  !print*, actbound, model_grid(I)%rho
!  bound = actbound
!  nOfPoints(I) = 0
! END DO
! ! now we will compute given numbers of points for the given spheres
! DO I = 1, Nvirtpart
!  point = ran2(idum)
!  DO J = 1, n_modelgrid
!   IF((bounds(J) > point)) THEN
!    nOfPoints(J) = nOfPoints(J) + 1
!    EXIT
!   END IF
!  END DO
! END DO
 ! now we will compute given numbers of points for the given spheres
 DO I = 1, n_modelgrid
  nOfPoints(I) = INT(FLOAT(Nvirtpart) * (model_grid(I)%rwind / R_inf) ** delta / sumr)
 END DO
 ! printing number of points for each model grid
 OPEN(UNIT=8,FILE='vp_distribution.dat')
  DO I = 1, n_modelgrid
   write(8,*) I, nOfPoints(I)
  END DO
 CLOSE(8)
 ! we have zero particles located
 NP = 0
 ! distribution of particles on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  np_shell = nOfPoints(I)
  !IF (np_shell == 0) STOP 'number of virtual particles is small'
  DO J = 1, np_shell
   NP = NP + 1
   CALL random_unitvector(direction)
   virtual_particle(NP)%pos = radius * direction
    write(20,*) virtual_particle(NP)%pos(1), virtual_particle(NP)%pos(2), virtual_particle(NP)%pos(3)
  END DO
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! 2D MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! for the case of 1D model
! we consider radial symmetric model
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 ALLOCATE (virtual_particle(Nvirtpart))
 print*, 'number of particles: ', Nvirtpart
 print*, 'computing positions of virtual particles...'
 rhomax = MAXVAL(model_grid(:)%rho)
 NP = 0
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 sumr = 0.D0
 DO I = 1, n_modelgrid
  sumr = sumr + (model_grid(I)%rwind / R_inf) ** delta
 END DO
 DO I = 1, n_modelgrid
  nOfPoints(I) = INT(FLOAT(Nvirtpart) * (model_grid(I)%rwind / R_inf) ** delta / sumr)
 END DO
 ! printing number of points for each model grid
 OPEN(UNIT=8,FILE='vp_distribution.dat')
  DO I = 1, n_modelgrid
   write(8,*) I, nOfPoints(I)
  END DO
 CLOSE(8)
 ! now we will compute given numbers of points for the given spheres
 DO I = 1, Nvirtpart
  point = 1.D2 * ran2(idum)
  DO J = 1, n_modelgrid
   IF((bounds(J) > point)) THEN
    nOfPoints(J) = nOfPoints(J) + 1
    EXIT
   END IF
  END DO
 END DO
 ! printing number of points for each model grid
! DO I = 1, n_modelgrid
!  print*, nOfPoints(I)
! END DO
! DO I = 1, n_modelgrid
!  nOfPoints(I) = 0
! END DO
! DO WHILE(NP < Npart)
!  rand_r = ran2(idum)
!  rand_rho = rhomax * ran2(idum + 1)
!  loc_rad = INT(FLOOR(1 + (n_modelgrid - 1) * rand_rho))
!!  print*, 'virtual_particles: rand_r = ', rand_r, ' rand_rho = ', rand_rho, ' loc_rad = ', loc_rad
!  IF(rand_rho <= model_grid(loc_rad)%rho) THEN
!   nOfPoints(loc_rad) = nOfPoints(loc_rad) + 1
!   NP = NP + 1
!   print*, 'found NP = ', NP
!  END IF
! END DO
! 
 ! we have zero particles located
 NP = 0
 ! distribution of particles on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = sqrt(model_grid(I)%rwind**2 - model_grid(I)%zwind**2)
  np_shell = nOfPoints(I)
  !IF (np_shell == 0) STOP 'number of virtual particles is small'
  DO J = 1, np_shell
   NP = NP + 1
   phi = 2.D0*pi*ran2(idum)
   virtual_particle(NP)%pos(1) = radius * cos(phi)
   virtual_particle(NP)%pos(2) = radius * sin(phi)
   virtual_particle(NP)%pos(3) = model_grid(I)%zwind
    write(20,*) virtual_particle(NP)%pos(1), virtual_particle(NP)%pos(2), virtual_particle(NP)%pos(3)
  END DO
 END DO

CASE DEFAULT
 STOP 'wrong choice of input model dimension...'
END SELECT
 CLOSE(20)

END SUBROUTINE
