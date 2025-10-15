! Scatter photon randomly from the photosphere isotropicaly in phi and
! with distribution function theta*d(theta) in theta
! INPUT: NONE
! OUTPUT: direction(DOUBLE(const_dimofspace)): calculated random direciton
! 
! 1x RETURN point
!
SUBROUTINE random_unitvector2(direction)
     
USE types
USE constants


IMPLICIT NONE 

DOUBLE PRECISION                 :: sint, cost, sinp, cosp, phi, ran2
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: direction

cost=ran2(idum)
cost=SQRT(cost)
sint=SQRT(1.D0 - cost*cost)
phi=2.D0*const_pi*ran2(idum)
cosp=COS(phi)
sinp=SIN(phi)
direction(ind_x)=sint*cosp
direction(ind_y)=sint*sinp                       
direction(ind_z)=cost 


! RETURN point
RETURN 

END SUBROUTINE random_unitvector2
