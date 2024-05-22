MODULE dummypacket

IMPLICIT NONE
SAVE

INTEGER, PARAMETER              :: n_dummy_packs = 10

TYPE dummyphoton 
   INTEGER                         :: cell_numb, pack_numb, active
   DOUBLE PRECISION                :: e_cmf, e_rf, freq_cmf, freq_rf, delta_s
   INTEGER                         :: typ, next_cross, last_line
   INTEGER                         :: n_interactions
   DOUBLE PRECISION, DIMENSION(3)  :: pos, dir 
   INTEGER                         :: l_ele, l_ion, l_lev, n_int = 0
   LOGICAL                         :: redShift, virtual
   LOGICAL                         :: occupied=.false.
END TYPE dummyphoton



TYPE(dummyphoton), DIMENSION(n_dummy_packs)             :: dummypackage



























END MODULE dummypacket
