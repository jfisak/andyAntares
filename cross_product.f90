! calculates cross product of two vectors
!
! INPUT: vec_a(DBLE(const_dimofspace)): the first vector
!        vec_b(DBLE(const_dimofspace)): the second vector
! OUTPUT: vysledek(DBLE(const_dimofspace)): the resulting vector
! 
SUBROUTINE cross_product(vec_a, vec_b, vysledek)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)                   :: vysledek
DOUBLE PRECISION, DIMENSION(const_dimofspace)                   :: vec_a, vec_b

vysledek(ind_x) = vec_a(ind_y) * vec_b(ind_z) - vec_a(ind_z) * vec_b(ind_y)
vysledek(ind_y) = vec_a(ind_z) * vec_b(ind_x) - vec_a(ind_x) * vec_b(ind_z)
vysledek(ind_z) = vec_a(ind_x) * vec_b(ind_y) - vec_a(ind_y) * vec_b(ind_x)


END SUBROUTINE
