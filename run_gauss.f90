  SUBROUTINE run_gauss(RG1)
!
! Two Gaussian random numbers generated from two uniform random
! numbers. Copyright (c) Tao Pang 1997.
!
  USE types

  IMPLICIT NONE
  
  DOUBLE PRECISION, INTENT (OUT) :: RG1
  DOUBLE PRECISION               :: R1, R2, ran2

  R1=-LOG(ran2(idum))
  R1=SQRT(R1)
  R2=2.D0*pi*ran2(idum)
  RG1=R1*SIN(R2)

  END SUBROUTINE run_gauss    
