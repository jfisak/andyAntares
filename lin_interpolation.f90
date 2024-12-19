SUBROUTINE lin_interpolation(pos0, vec1, pos1, vec2, pos2, int_vector)
! SUBROUTINE lin_interpolation(pos0, vec1, pos1, vec2, pos2, int_vector, vecpos0, neighborscells)

USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)                 :: vec1, vec2
DOUBLE PRECISION                               :: pos0, pos1, pos2
DOUBLE PRECISION, DIMENSION(const_dimofspace)                 :: int_vector
DOUBLE PRECISION, DIMENSION(const_dimofspace)                 :: lina, linb
! DOUBLE PRECISION, DIMENSION(3)                 :: vecpos0
! INTEGER, DIMENSION(8)                          :: neighborscells

if(pos1 == pos2) then
 write(*,*) 'lin_interpolation: pos1 = ', pos1, ' pos2 = ', pos2
 STOP 'pos1 == pos2'
end if

IF(pos1 > pos2 .and. (pos0 < pos2 .or. pos0 > pos1)) then
 write(*,*) 'lin_interpolation: pos1 = ', pos1, ' pos0 = ', pos0, ' pos2 = ', pos2
 ! DO I = 1,8
 !  write(3,*) dyn_cell(neighborscells(I))%corner, dyn_cell(neighborscells(I))%width
 ! END DO
 ! write(4,*) vecpos0
 STOP 'pos0 is not in the interval'
ELSE IF(pos1 < pos2 .and. (pos0 > pos2 .or. pos0 < pos1)) THEN
 write(*,*) 'lin_interpolation: pos1 = ', pos1, ' pos0 = ', pos0, ' pos2 = ', pos2
 ! DO I = 1,8
 !  write(3,*) dyn_cell(neighborscells(I))%corner, dyn_cell(neighborscells(I))%width
 ! END DO
 ! write(4,*) vecpos0
 STOP 'pos0 is not in the interval'
END IF

lina = (vec2-vec1)/(pos2-pos1)
linb = -(vec2*pos1-vec1*pos2)/(pos2-pos1)
int_vector = lina*pos0 + linb

! write(*,*) 'lin_interpolation: vec1 = ', vec1, ' vec2 = ', vec2
! write(*,*) 'lin_interpolation: pos2 = ', pos2, ' pos1 = ', pos1
! write(*,*) 'lin_interpolation: vec1 = ', norm2(vec1)/const_c, ' vec2 = ', norm2(vec2)/const_c
! write(*,*) 'lin_interpolation: int_vector = ', int_vector
! write(*,*) 'lin_interpolation: int_vector = ', norm2(int_vector)/const_c

END SUBROUTINE lin_interpolation
