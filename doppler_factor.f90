SUBROUTINE doppler_factor(pack_index, D)

! Calculate the Doppler factor to transform a rest frame frequency
! into a comoving frame frequency (Mihalas & Mihalas Eq. 89.5)

USE types

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: D_gamma, D
DOUBLE PRECISION, DIMENSION(3)    :: vel_vec
  
D_gamma = 1.D0 ! For non-relativistic case    
 CALL velo(pack_index, vel_vec, velApprox)
D = D_gamma * (1.D0 -  DOT_PRODUCT(package(pack_index)%dir,vel_vec)/light_speed)

END SUBROUTINE doppler_factor
