! this sbr will be activated when the plasma in the current propagation cell is optically thick and
! Monte Carlo is not efficient instead of difficult propagation a new packet is radiated at the
! propGrid cell boundary
SUBROUTINE do_dpackage(pack_index)

USE types
IMPLICIT NONE

INTEGER                                 :: pack_index
INTEGER                                 :: cur_pgi


cur_pgi = package(pack_index)%cell_numb









END SUBROUTINE do_dpackage
