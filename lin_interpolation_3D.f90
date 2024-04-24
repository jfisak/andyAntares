SUBROUTINE lin_interpolation_3D(pos_0, pos_1, pos_2, vel_1, vel_2, cur_coo, vec_vel)

USE TYPES
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                                  :: pos_0, pos_1, pos_2, vel_1, vel_2
DOUBLE PRECISION, DIMENSION(3)                                  :: vec_vel
DOUBLE PRECISION, DIMENSION(3)                                  :: point_pos
INTEGER                                                         :: cur_coo
DOUBLE PRECISION, DIMENSION(3)                                  :: lina, linb


point_pos = pos_1 + (pos_2 - pos_1) * (pos_0(cur_coo) - pos_1(cur_coo))/(pos_2(cur_coo) - pos_1(cur_coo))

! the interpolation must be processed for each dimension separately
lina = (pos_0(cur_coo) - pos_1(cur_coo))/(pos_1(cur_coo)-pos_2(cur_coo)) * (vel_2 - vel_1)
linb = vel_1

vec_vel = lina + linb

IF(isnan(vec_vel(1)) .or. isnan(vec_vel(2)) .or. isnan(vec_vel(3))) THEN
 write(*,*) 'lin_interpolation_3D: pos_1 = ', pos_1, ' pos_2 = ', pos_2
 write(*,*) 'lin_interpolation_3D: cur_coo = ', cur_coo
 write(*,*) 'lin_interpolation_3D: point_pos = ', point_pos, ' lina = ', lina, ' linb = ', linb
 STOP 'lin_interpolation_3D the velocity is NaN'
END IF


END SUBROUTINE lin_interpolation_3D
