SUBROUTINE vel_discrete_points(pack_index, vel_vec)

IMPLICIT NONE

USE types

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec

! spherically symmetric model
IF(model_type == 1) THEN
 

END IF



END SUBROUTINE vel_discrete_points
