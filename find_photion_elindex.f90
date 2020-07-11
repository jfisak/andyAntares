! this SBR will find for the given photoionisation cross section data
! the corresponding atomic energy level
SUBROUTINE find_photion_elindex(indexe, indexi, index_level, phot_index)

USE types

IMPLICIT NONE

INTEGER                                 :: indexe, indexi
INTEGER                                 :: index_level
INTEGER                                 :: phot_index
INTEGER                                 :: nlevs
DOUBLE PRECISION, PARAMETER             :: rydberg = 13.5979996 !(eV)
INTEGER                                 :: indexl
INTEGER                                 :: cur_indexl

nlevs = SIZE(elements(indexe)%ions(indexi)%levels)

phot_index = -1
! write(*,*) 'find_photion_elindex: index_level = ', index_level
DO indexl = 1, nlevs
 cur_indexl = elements(indexe)%ions(indexi)%levels(indexl)%levelindex
 IF(index_level == cur_indexl) THEN
  phot_index = indexl
  EXIT
 END IF
END DO
! write(*,*) 'find_photion_elindex: phot_index = ', phot_index



END SUBROUTINE
