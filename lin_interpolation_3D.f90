! calculates an interpolation between two points in full 3D space
!
! INPUT: pos_0(DBLE(const_dimofspace)): position where the interpolation is processed
!        pos_1(DBLE(const_dimofspace)): the first interpolation point
!        pos_2(DBLE(const_dimofspace)): the second interpolation point
!        vel_1(DBLE(const_dimofspace)): velocity vector in the first interpolation point
!        vel_2(DBLE(const_dimofspace)): velocity vector in the second interpolation point
!        cur_coo(INT): index of the coordinate where the interpolation is processed
! OUTPUT: vec_vel(DBLE(const_dimofspace)): calculated velocity vector
!
SUBROUTINE lin_interpolation_3D(pos_0, pos_1, pos_2, vel_1, vel_2, cur_coo, vec_vel)

USE TYPES
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)                           :: pos_0, pos_1, pos_2, vel_1, vel_2
DOUBLE PRECISION, DIMENSION(const_dimofspace)                           :: vec_vel
DOUBLE PRECISION, DIMENSION(const_dimofspace)                           :: point_pos
INTEGER                                                                 :: cur_coo
DOUBLE PRECISION, DIMENSION(const_dimofspace)                           :: lina, linb


point_pos = pos_1 + (pos_2 - pos_1) * (pos_0(cur_coo) - pos_1(cur_coo))/(pos_2(cur_coo) - pos_1(cur_coo))

! the interpolation must be processed for each dimension separately
lina = (pos_0(cur_coo) - pos_1(cur_coo))/(pos_1(cur_coo)-pos_2(cur_coo)) * (vel_2 - vel_1)
linb = vel_1

vec_vel = lina + linb

IF(isnan(vec_vel(ind_x)) .or. isnan(vec_vel(ind_y)) .or. isnan(vec_vel(ind_z))) THEN
 write(*,*) 'lin_interpolation_3D: pos_1 = ', pos_1, ' pos_2 = ', pos_2
 write(*,*) 'lin_interpolation_3D: cur_coo = ', cur_coo
 write(*,*) 'lin_interpolation_3D: point_pos = ', point_pos, ' lina = ', lina, ' linb = ', linb
 STOP 'lin_interpolation_3D the velocity is NaN'
END IF


END SUBROUTINE lin_interpolation_3D
