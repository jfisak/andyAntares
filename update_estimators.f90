SUBROUTINE update_estimators(pack_index, dist)

  USE types

  IMPLICIT NONE

    INTEGER                           :: pack_index, current_mgi, get_package_model_index
    DOUBLE PRECISION                  :: dist
 
    current_mgi = get_package_model_index(pack_index)

    model_grid(current_mgi)%J = model_grid(current_mgi)%J + package(pack_index)%e_cmf * dist

END SUBROUTINE update_estimators
