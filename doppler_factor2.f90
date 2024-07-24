SUBROUTINE doppler_factor2(pack_index, cur_pos, cur_mgi, D)

! Calculate the Doppler factor to transform a rest frame frequency
! into a comoving frame frequency (Mihalas & Mihalas Eq. 89.5)

USE types
USE constants
USE dummypacket

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: D_gamma, D
DOUBLE PRECISION, DIMENSION(const_dimofspace)    :: vel_vec

DOUBLE PRECISION, DIMENSION(const_dimofspace)    :: cur_pos, cur_dir
INTEGER                                 :: cur_mgi, cur_dummy_index
  
IF(pack_index < SIZE(package)) THEN
 cur_dir = package(pack_index)%dir
ELSE IF(pack_index > SIZE(package)) THEN
 cur_dummy_index = pack_index - SIZE(package)
 cur_dir = dummypackage(cur_dummy_index)%dir
END IF
D_gamma = 1.D0 ! For non-relativistic case    
CALL velo(pack_index, cur_mgi, vel_vec, velApprox)
D_gamma = 1/sqrt(1-norm2(vel_vec)**2.0/light_speed**2.00)
D = D_gamma * (1.D0 -  DOT_PRODUCT(cur_dir,vel_vec)/light_speed)

END SUBROUTINE doppler_factor2
