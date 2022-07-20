SUBROUTINE lin_interpolation(pos0, vec1, pos1, vec2, pos2, int_vector, vecpos0, neighborscells)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                 :: vec1, vec2
DOUBLE PRECISION                               :: pos0, pos1, pos2
DOUBLE PRECISION, DIMENSION(3)                 :: int_vector
DOUBLE PRECISION, DIMENSION(3)                 :: lina, linb
DOUBLE PRECISION, DIMENSION(3)                 :: vecpos0
INTEGER, DIMENSION(8)                          :: neighborscells
INTEGER                                        :: pack_index, I

IF(pos1 > pos2 .and. (pos0 < pos2 .or. pos0 > pos1)) then
 write(*,*) 'lin_interpolation: pack_index = ', pack_index
 write(*,*) 'lin_interpolation: pos0 = ', pos0, ' pos1 = ', pos1, ' pos2 = ', pos2
 write(*,*) 'lin_interpolation: pos0 = ', pos0 - pos1, ' pos2 = ', pos0 - pos2
 DO I = 1,8
  write(3,*) dyn_cell(neighborscells(I))%corner, dyn_cell(neighborscells(I))%width
 END DO
 write(4,*) vecpos0
 STOP 'pos0 is not in the interval'
ELSE IF(pos1 < pos2 .and. (pos0 > pos2 .or. pos0 < pos1)) THEN
 write(*,*) 'lin_interpolation: pack_index = ', pack_index
 write(*,*) 'lin_interpolation: pos0 = ', pos0, ' pos1 = ', pos1, ' pos2 = ', pos2
 write(*,*) 'lin_interpolation: pos0 = ', pos0 - pos1, ' pos2 = ', pos0 - pos2
 DO I = 1,8
  write(3,*) dyn_cell(neighborscells(I))%corner, dyn_cell(neighborscells(I))%width
 END DO
 write(4,*) vecpos0
 STOP 'pos0 is not in the interval'
END IF

lina = (vec2-vec1)/(pos2-pos1)
linb = -(vec2*pos1-vec1*pos2)/(pos2-pos1)
int_vector = lina*pos0 + linb

! write(*,*) 'lin_interpolation: vec1 = ', vec1, ' vec2 = ', vec2
! write(*,*) 'lin_interpolation: vec1 = ', norm2(vec1)/light_speed, ' vec2 = ', norm2(vec2)/light_speed
! write(*,*) 'lin_interpolation: int_vector = ', int_vector
! write(*,*) 'lin_interpolation: int_vector = ', norm2(int_vector)/light_speed

END SUBROUTINE lin_interpolation
