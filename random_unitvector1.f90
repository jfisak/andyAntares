! calculates random unitvector in an arbitrary direction
!
! INPUT: sint(DBL): sin(theta)
!        cost(DBL): cos(theta)
!        sinp(DBL): sin(phi)
!        cosp(DBL): cos(phi)
! OUTPUT: direction(DBL(const_dimofspace)): random direction
!
! 1x RETURN point
!
  SUBROUTINE random_unitvector1(direction, sint, cost, sinp, cosp)

!     Scatter photon randomly and isotropicaly (in phi and theta) in any direction in the wind
      USE types
USE constants

      IMPLICIT NONE 

      DOUBLE PRECISION                 :: sint, cost, sinp, cosp, phi, ran2
      DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: direction
!      INTEGER                          :: idum
!      COMMON / RAN_SEED / idum

      cost=2.D0 * ran2(idum) - 1.D0
      sint=SQRT(1.D0 - cost*cost)
      phi=2.D0*const_pi*ran2(idum)
      cosp=COS(phi)
      sinp=SIN(phi)
      direction(ind_x)=sint*cosp
      direction(ind_y)=sint*sinp                       
      direction(ind_z)=cost  

!      length=SQRT(direction(1)**2 + direction(2)**2 + direction(3)**2)  
!      print*,length
!      direction(1)=direction(1)/length
!      direction(2)=direction(2)/length
!      direction(3)=direction(3)/length

!      dir=SQRT(norm_direction(1)**2 + norm_direction(2)**2 + norm_direction(3)**2)  
!      print*, dir
!      STOP
      ! RETURN point
      RETURN 

  END SUBROUTINE random_unitvector1
