!  Scatter photon randomly and isotropicaly (in phi and theta) in any direction in the wind
!
! INPUT: NONE
! OUTPUT: direction(DBLE(const_dimofspace)): a calculated random direction
!
! 1x RETURN point
!
SUBROUTINE random_unitvector(direction)

USE types
USE constants

IMPLICIT NONE 

DOUBLE PRECISION                 :: sint, cost, sinp, cosp, phi, ran2
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: direction
!      INTEGER                          :: idum
!      COMMON / RAN_SEED / idum

cost=2.D0*ran2(idum) - 1.D0
sint=SQRT(1.D0 - cost*cost)
!phi=2.D0*const_pi*ran2(idum)
phi=2.D0*const_pi*ran2(idum)
cosp=COS(phi)
sinp=SIN(phi)
direction(ind_x)=sint*cosp
direction(ind_y)=sint*sinp                       
direction(ind_z)=cost  

! RETURN point
RETURN 

END SUBROUTINE random_unitvector
