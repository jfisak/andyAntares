
 SUBROUTINE velo(pack_index,vel_vec)

  USE types

  IMPLICIT NONE    

    INTEGER                           :: pack_index
    DOUBLE PRECISION                  :: vel_radial, vec_length
    DOUBLE PRECISION, DIMENSION(3)    :: vel_vec

!    print*, R_inf, V_inf
!    vel_radial = V_inf/R_inf * vec_length(package(pack_index)%pos)
   vel_radial = V_inf * (1.D0 - b/vec_length(package(pack_index)%pos))**beta

    vel_vec = package(pack_index)%pos/vec_length(package(pack_index)%pos) * vel_radial

    !print*, vec_length(package(pack_index)%pos), vel_radial, vec_length(vel_vec), V_inf, R_inf/r_sun

 END SUBROUTINE velo
