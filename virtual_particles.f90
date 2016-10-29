! this subroutine creates virtual particles for a dynamical grid calculation
SUBROUTINE virtual_particles(dimIM)

 USE types

 IMPLICIT NONE 

 ! input variables
 INTEGER                        :: dimIM
 ! 
 INTEGER                        :: I,J,K,NP
 ! number of particles
 INTEGER, PARAMETER             :: Npart=100000
 ! 1D model: intervals for particles distribution
 INTEGER                        :: Ntheta, Nphi
 DOUBLE PRECISION               :: radius, phi, theta
 DOUBLE PRECISION               :: angle
 DOUBLE PRECISION, DIMENSION(3) :: direction
 INTEGER                        :: angleParam, np_shell
 ! division of an interval [0, 1] into parts corresponding to a density
 DOUBLE PRECISION               :: rhotot
 ! bound of the division
 DOUBLE PRECISION               :: bound, actbound
 DOUBLE PRECISION, DIMENSION(n_modelgrid) :: bounds
 INTEGER, DIMENSION(n_modelgrid) :: nOfPoints
 ! a random point
 DOUBLE PRECISION               :: point
 DOUBLE PRECISION               :: ran2
 ! not needed variables, only for a subroutine call
 DOUBLE PRECISION               :: ct, st, cp, sp


 SELECT CASE (dimIM)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!! 1D MODEL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! for the case of 1D model
 ! we consider radial symmetric model
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE(1)
 OPEN(20,FILE="virtual_particles.dat")
 ! angle between particles
 ! angleParam = 20
 ! number of virtual particles
 ! Ntheta =  angleParam
 ! Nphi = 2* angleParam
 !Npart = n_modelgrid * (Ntheta - 2) * Nphi + 2 * n_modelgrid
 print*, 'number of particles: ', Npart
 print*, 'computing positions of virtual particles...'
 ALLOCATE (virtual_particle(Npart))
 ! computing number of points on a shell from a density
 ! firstly we compute a total number of density
 rhotot = 0.D0
 DO I = 1, n_modelgrid
 rhotot = rhotot + model_grid(I)%rho
 END DO
 ! now we divide an interval [0, 1] into parts which lenght
 ! corresponds to the density magnitude
 bound = 0.D0
 DO I = 1, n_modelgrid
  actbound = bound + model_grid(I)%rho / rhotot
  bounds(I) = actbound
  !print*, actbound, model_grid(I)%rho
  bound = actbound
  nOfPoints(I) = 0
 END DO
 ! now we will compute given numbers of points for the given spheres
 DO I = 1, Npart
  point = ran2(idum)
  DO J = 1, n_modelgrid
   IF((bounds(J) > point)) THEN
    nOfPoints(J) = nOfPoints(J) + 1
    EXIT
   END IF
  END DO
 END DO
 ! printing number of points for each model grid
 DO I = 1, n_modelgrid
  print*, nOfPoints(I)
 END DO
 ! we have zero particles located
 NP = 0
 ! distribution of particles on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  np_shell = nOfPoints(I)
  IF (np_shell == 0) STOP 'number of virtual particles is small'
  DO J = 1, np_shell
   NP = NP + 1
   CALL random_unitvector1(direction,st,ct,sp,cp)
   virtual_particle(NP)%pos = radius * direction
    write(20,*) virtual_particle(NP)%pos(1), virtual_particle(NP)%pos(2), virtual_particle(NP)%pos(3)
  END DO
 END DO
 CLOSE(20)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 2D models
 CASE(2)
 CASE DEFAULT
  STOP 'wrong choice of input model dimension...'
 END SELECT

END SUBROUTINE
