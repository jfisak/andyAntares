SUBROUTINE lin_interpolation(pos0, vec1, pos1, vec2, pos2, int_vector)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                 :: vec1, vec2
DOUBLE PRECISION                               :: pos0, pos1, pos2
DOUBLE PRECISION, DIMENSION(3)                 :: int_vector
DOUBLE PRECISION, DIMENSION(3)                 :: lina, linb

lina = (vec2-vec1)/(pos2-pos1)
linb = (vec2*pos1-vec1*pos2)/(pos2-pos1)
int_vector = lina*pos0 + linb

END SUBROUTINE lin_interpolation
