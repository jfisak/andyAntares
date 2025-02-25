! calculates Monte Carlo estimators
!

SUBROUTINE update_estimators(pack_index, dist)

USE types
USE constants

IMPLICIT NONE

INTEGER                           :: pack_index, current_mgi, get_package_model_index
DOUBLE PRECISION                  :: dist
DOUBLE PRECISION, PARAMETER       :: delta = 1.D-1
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_pos, cur_dir, rad_dir
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: deltas

current_mgi = get_package_model_index(pack_index)

model_grid(current_mgi)%J = model_grid(current_mgi)%J + package(pack_index)%e_cmf * dist


! I in radial direction
cur_pos = package(pack_index)%pos
cur_dir = package(pack_index)%dir
rad_dir = cur_pos/norm2(cur_pos)
deltas = abs(cur_dir(:) - rad_dir(:))
IF(maxval(deltas(:)) < delta) THEN
 model_grid(current_mgi)%Irad = model_grid(current_mgi)%Irad + package(pack_index)%e_cmf * dist
END IF



model_grid(current_mgi)%Frad = model_grid(current_mgi)%Frad + &
 cur_dir * package(pack_index)%e_cmf * dist



END SUBROUTINE update_estimators
