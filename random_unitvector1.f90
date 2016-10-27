  SUBROUTINE random_unitvector1(direction, sint, cost, sinp, cosp)

!     Scatter photon randomly and isotropicaly (in phi and theta) in any direction in the wind
      USE types

      IMPLICIT NONE 

!      DOUBLE PRECISION, PARAMETER      :: pi=3.1415926535897932D0 
      DOUBLE PRECISION                 :: sint, cost, sinp, cosp, phi, ran2, length
      DOUBLE PRECISION, DIMENSION(3)   :: direction
!      INTEGER                          :: idum
!      COMMON / RAN_SEED / idum

      cost=2.D0*ran2(idum) - 1.D0
      sint=SQRT(1.D0 - cost*cost)
! it seems not to be correct because this distribution is not uniform -- points are
! more distributed arround the poles
!      phi=2.D0*pi*ran2(idum)
      phi = acos(2.0*ran2(idum)-1.0)
      cosp=COS(phi)
      sinp=SIN(phi)
      direction(1)=sint*cosp
      direction(2)=sint*sinp                       
      direction(3)=cost  

!      length=SQRT(direction(1)**2 + direction(2)**2 + direction(3)**2)  
!      print*,length
!      direction(1)=direction(1)/length
!      direction(2)=direction(2)/length
!      direction(3)=direction(3)/length

!      dir=SQRT(norm_direction(1)**2 + norm_direction(2)**2 + norm_direction(3)**2)  
!      print*, dir
!      STOP
      RETURN 

  END SUBROUTINE random_unitvector1
