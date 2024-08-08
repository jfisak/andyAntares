MODULE dummypacket

IMPLICIT NONE
SAVE

INTEGER, PARAMETER              :: n_dummy_packs = 10

TYPE dummyphoton 
   INTEGER                         :: cell_numb, active
   DOUBLE PRECISION                :: e_cmf, e_rf, freq_cmf, freq_rf, delta_s
   INTEGER                         :: typ, next_cross, last_line
   INTEGER                         :: n_interactions
   DOUBLE PRECISION, DIMENSION(3)  :: pos, dir 
   INTEGER                         :: l_ele, l_ion, l_lev, n_int = 0
   LOGICAL                         :: redShift, virtual
   LOGICAL                         :: occupied=.false.
END TYPE dummyphoton


 
TYPE(dummyphoton), DIMENSION(n_dummy_packs)             :: dummypackage

CONTAINS

 ! a function that seeks the free index
 INTEGER function find_free_index()
  IMPLICIT NONE
  INTEGER                        :: cur_index
  LOGICAL                        :: is_occup
 
  DO cur_index = 1, n_dummy_packs
   is_occup = dummypackage(cur_index)%occupied
   IF(.NOT. is_occup) THEN
    dummypackage(cur_index)%occupied = .true.
    find_free_index = cur_index
    EXIT
   END IF
  END DO
  IF(is_occup .and. cur_index == n_dummy_packs) THEN
   write(*,*) 'dummypacket, find_free_index: no free index found!'
   write(*,*) 'probably some dummy packages have not been deactivated'
   STOP 'dummypacket, find_free_index'
  END IF
 end function

 SUBROUTINE deactivate_dummy_packet(pack_index)
  USE types

  IMPLICIT NONE

  INTEGER                               :: pack_index

  dummypackage(pack_index)%occupied = .false.
  dummypackage(pack_index)%pos = (/ 0.D0, 0.D0, 0.D0 /)
  dummypackage(pack_index)%dir = (/ 0.D0, 0.D0, 0.D0 /)
  dummypackage(pack_index)%cell_numb = 0
  dummypackage(pack_index)%freq_cmf = 0.D0
  dummypackage(pack_index)%freq_rf = 0.D0
  dummypackage(pack_index)%e_cmf = 0.D0
  dummypackage(pack_index)%e_rf = 0.D0
  dummypackage(pack_index)%typ = NONE
  dummypackage(pack_index)%next_cross = 0
  dummypackage(pack_index)%last_line = NONE


 END SUBROUTINE

 ! copy a normal package to a dummy packet
 SUBROUTINE copy_package(input_packet, ind_dummy_packet)
  USE types

  IMPLICIT NONE

  INTEGER                               :: input_packet, ind_dummy_packet

  dummypackage(ind_dummy_packet)%pos = package(input_packet)%pos
  dummypackage(ind_dummy_packet)%dir = package(input_packet)%dir
  dummypackage(ind_dummy_packet)%cell_numb = package(input_packet)%cell_numb
  dummypackage(ind_dummy_packet)%freq_cmf = package(input_packet)%freq_cmf
  dummypackage(ind_dummy_packet)%freq_rf = package(input_packet)%freq_rf
  dummypackage(ind_dummy_packet)%e_cmf = package(input_packet)%e_cmf
  dummypackage(ind_dummy_packet)%e_rf = package(input_packet)%e_rf
  dummypackage(ind_dummy_packet)%typ = package(input_packet)%typ
  dummypackage(ind_dummy_packet)%next_cross = package(input_packet)%next_cross
  dummypackage(ind_dummy_packet)%last_line = package(input_packet)%last_line

 END SUBROUTINE

 SUBROUTINE teleport_dummypacket(pack_index, new_pos)
  USE types

  IMPLICIT NONE

  INTEGER                                               :: pack_index
  DOUBLE PRECISION, DIMENSION(const_dimofspace)         :: new_pos


  ! setting a new position
  dummypackage(pack_index)%pos = new_pos

 END SUBROUTINE

END MODULE dummypacket
