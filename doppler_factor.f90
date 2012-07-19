 SUBROUTINE doppler_factor(pack_index, D)
 
  USE types

  IMPLICIT NONE    

    INTEGER                           :: pack_index
    DOUBLE PRECISION                  :: D_gamma, D, vel_radial, vec_length
    DOUBLE PRECISION, DIMENSION(3)    :: vel_vec
    
    D_gamma = 1.D0 ! For non-relativistic case    
    CALL velo(pack_index, vel_vec)
    D = D_gamma * (1.D0 -  DOT_PRODUCT(package(pack_index)%dir,vel_vec)/light_speed)

 END SUBROUTINE doppler_factor
