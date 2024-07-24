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
    find_free_index = cur_index
    EXIT
   END IF
  END DO
 end function

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



END MODULE dummypacket
