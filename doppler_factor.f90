SUBROUTINE doppler_factor(pack_index, D)

! Calculate the Doppler factor to transform a rest frame frequency
! into a comoving frame frequency (Mihalas & Mihalas Eq. 89.5)

USE types
USE constants
USE dummypacket

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: D_gamma, D
DOUBLE PRECISION, DIMENSION(const_dimofspace)    :: vel_vec, cur_dir

INTEGER                                 :: cur_mgi, get_package_model_index
INTEGER                                 :: dummypack_index
  
cur_mgi = get_package_model_index(pack_index)
IF(pack_index <= SIZE(package)) THEN
 cur_dir = package(pack_index)%dir
 CALL velo(pack_index, vel_vec, 0)
ELSE IF(pack_index > SIZE(package)) THEN
 dummypack_index = pack_index - SIZE(package)
 cur_dir = dummypackage(dummypack_index)%dir
 CALL velo(dummypack_index, vel_vec, 0)
END IF


D_gamma = 1.D0 ! For non-relativistic case    
! for the relativistic case
! D_gamma = 1/sqrt(1-norm2(vel_vec)**2.0/const_c**2.00)
D = D_gamma * (1.D0 -  DOT_PRODUCT(cur_dir,vel_vec)/const_c)


END SUBROUTINE doppler_factor
