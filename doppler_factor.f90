SUBROUTINE doppler_factor(pack_index, D)

! Calculate the Doppler factor to transform a rest frame frequency
! into a comoving frame frequency (Mihalas & Mihalas Eq. 89.5)

USE types
USE constants

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: D_gamma, D
DOUBLE PRECISION, DIMENSION(3)    :: vel_vec

DOUBLE PRECISION, DIMENSION(3)          :: cur_pos
INTEGER                                 :: cur_mgi, get_package_model_index
  
cur_pos = package(pack_index)%pos
cur_mgi = get_package_model_index(pack_index)

D_gamma = 1.D0 ! For non-relativistic case    
 CALL velo(pack_index, cur_pos, cur_mgi, vel_vec, velApprox)
D_gamma = 1/sqrt(1-norm2(vel_vec)**2.0/light_speed**2.00)
D = D_gamma * (1.D0 -  DOT_PRODUCT(package(pack_index)%dir,vel_vec)/light_speed)


END SUBROUTINE doppler_factor
