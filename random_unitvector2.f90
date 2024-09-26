  SUBROUTINE random_unitvector2(direction)
     
      USE types
USE constants

!     Scatter photon randomly from the photosphere isotropicaly in phi and with distribution function theta*d(theta) in theta

      IMPLICIT NONE 

      DOUBLE PRECISION                 :: sint, cost, sinp, cosp, phi, ran2
      DOUBLE PRECISION, DIMENSION(3)   :: direction

      cost=ran2(idum)
      cost=SQRT(cost)
      sint=SQRT(1.D0 - cost*cost)
      phi=2.D0*const_pi*ran2(idum)
      cosp=COS(phi)
      sinp=SIN(phi)
      direction(1)=sint*cosp
      direction(2)=sint*sinp                       
      direction(3)=cost 


      RETURN 

  END SUBROUTINE random_unitvector2
