SUBROUTINE vel_discrete_points_constant(pack_index, vel_vec)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: vel_vec, pos_0

INTEGER                                 :: cur_mgi, cur_pgi
INTEGER                                 :: get_package_model_index

cur_mgi = get_package_model_index(pack_index)
! pos_0 = package(pack_index)%pos
cur_pgi = package(pack_index)%cell_numb

pos_0 = dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width / 2.0

IF(model_type == 1) THEN
 vel_vec = model_grid(cur_mgi)%vel * pos_0 / norm2(pos_0)
 write(*,*) 'vel_discrete_points: pack_index = ', pack_index, ' vel_vec = ', vel_vec
ELSE IF(model_type == 2) THEN
 vel_vec = model_grid(cur_mgi)%velocity
END IF


END SUBROUTINE vel_discrete_points_constant
