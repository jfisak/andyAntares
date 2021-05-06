SUBROUTINE vel_discrete_points_interpolation(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

INTEGER                                 :: get_package_model_index

cur_mgi = get_package_model_index(pack_index)

IF(model_type == 1) THEN
 vel_vec = model_grid(cur_mgi)%vel * pos_0 / norm2(pos_0)
 write(*,*) 'vel_discrete_points: vel_0 = ', vel_0, ' pos_0 = ', pos_0
ELSE IF(model_type == 2) THEN
 vel_vec = model_grid(cur_mgi)%velocity
END IF


END SUBROUTINE vel_discrete_points_interpolation
