SUBROUTINE emit_rpackage(pack_index)

  USE types

  IMPLICIT NONE    

    INTEGER                           :: pack_index
    DOUBLE PRECISION, DIMENSION(3)    :: direction

    package(pack_index)%typ = type_rpkt
    CALL random_unitvector(direction)
    package(pack_index)%dir = direction
    package(pack_index)%last_cross = NONE

! Frequency transformation from cmf to rf should be done later
    
END

