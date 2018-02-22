  SUBROUTINE random_unitvector2(direction)
     
      USE types

!     Scatter photon randomly from the photosphere isotropicaly in phi and with distribution function theta*d(theta) in theta

      IMPLICIT NONE 

!      DOUBLE PRECISION, PARAMETER      :: pi = 3.1415926535897932D0
      DOUBLE PRECISION                 :: sint, cost, sinp, cosp, phi, ran2, length
      DOUBLE PRECISION, DIMENSION(3)   :: direction
      REAL(8)                           :: random
!      INTEGER                          :: idum
!      COMMON / RAN_SEED / idum

      cost=DBLE(random())
      cost=SQRT(cost)
      sint=SQRT(1.D0 - cost*cost)
      phi=2.D0*pi*DBLE(random())
      cosp=COS(phi)
      sinp=SIN(phi)
      direction(1)=sint*cosp
      direction(2)=sint*sinp                       
      direction(3)=cost 

!      length=SQRT(direction(1)**2 + direction(2)**2 + direction(3)**2)  
!      direction(1)=direction(1)/length
!      direction(2)=direction(2)/length
!      direction(3)=direction(3)/length      

      RETURN 

  END SUBROUTINE random_unitvector2
