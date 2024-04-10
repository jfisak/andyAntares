SUBROUTINE lin_interpolation_3D(pos_0, pos_1, pos_2, vel_1, vel_2, cur_coo, vec_vel)

USE TYPES
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                                  :: pos_0, pos_1, pos_2, vel_1, vel_2
DOUBLE PRECISION, DIMENSION(3)                                  :: vec_vel
DOUBLE PRECISION, DIMENSION(3)                                  :: point_pos
INTEGER                                                         :: cur_coo




point_pos = pos_1 + (pos_2 - pos_1) * &
  & (pos_0(cur_coo) - pos_1(cur_coo))/(pos_1(cur_coo) - pos_2(cur_coo))







END SUBROUTINE lin_interpolation_3D
