! this subroutine creates virtual particles for a dynamical grid calculation
SUBROUTINE virtual_particles(dimIM)

 USE types

 IMPLICIT NONE 

 ! input variables
 INTEGER                        :: dimIM
 ! 
 INTEGER                        :: I,J,K,NP
 ! number of particles
 INTEGER, PARAMETER             :: Npart=50000
 ! 1D model: intervals for particles distribution
 INTEGER                        :: Ntheta, Nphi
 DOUBLE PRECISION               :: radius, phi, theta
 DOUBLE PRECISION               :: angle
 DOUBLE PRECISION, DIMENSION(3) :: direction
 INTEGER                        :: angleParam, np_shell


 SELECT CASE (dimIM)
 ! for the case of 1D model
 ! we consider radial symmetric model
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
 ! we have zero particles located
 NP = 0
 ! distribution of particles on the shell of the radius R
 DO I = 1, n_modelgrid
  radius = model_grid(I)%rwind
  np_shell = Npart/n_modelgrid
  IF (np_shell == 0) STOP 'number of virtual particles is small'
  DO J = 1, np_shell
   NP = NP + 1
   CALL random_unitvector2(direction)
   IF(mod(NP,2).EQ.0) direction(3) = -direction(3)
   virtual_particle(NP)%pos = radius * direction
    write(20,*) virtual_particle(NP)%pos(1), virtual_particle(NP)%pos(2), virtual_particle(NP)%pos(3)
  END DO
 END DO
 !DO I = 1, n_modelgrid
 ! radius = model_grid(I)%rwind
 ! ! we define number of particles in shells
 ! DO J = 1, Ntheta
 !  ! if the particles are in the pole we have to treat them separately
 !  IF( (J.EQ.1).OR.(J.EQ.Ntheta)) THEN
 !   virtual_particle(J)%pos(1) = 0.D0
 !   virtual_particle(J)%pos(2) = 0.D0
 !   if (J.EQ.1) virtual_particle(J)%pos(3) = -radius
 !   if (J.EQ.Ntheta) virtual_particle(J)%pos(3) = radius
 !   NP = NP + 1
 !   ! we do not have to calculate the phi positions of particles
 !   CYCLE
 !  END IF
 !  theta = (2.E0 * pi) / Ntheta * REAL(J) - pi
 !  DO K = 1, Nphi
 !   phi = (2 * pi)/(Nphi) * K
 !   NP = NP + 1
 !   virtual_particle(NP)%pos(1) = radius * cos(theta) * cos(phi)
 !   virtual_particle(NP)%pos(2) = radius * cos(theta) * sin(phi)
 !   virtual_particle(NP)%pos(3) = radius * sin(theta)
 !   write(20,*) virtual_particle(NP)%pos(1), virtual_particle(NP)%pos(2), virtual_particle(NP)%pos(3)
 !  END DO
 ! END DO
 !END DO
 CLOSE(20)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 2D models
 CASE(2)
 CASE DEFAULT
  STOP 'wrong choice of input model dimension...'
 END SELECT

END SUBROUTINE
