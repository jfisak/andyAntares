! corrects the position of a packet if its position is slightly different from the current
! propGrid cell index
!
! INPUT: cur_index(INT): x, y, z index
!        pack_index(INT): curent packet index
!        cur_pgi(INT): propGrid index
!
SUBROUTINE correction_propagation(cur_index, pack_index, cur_pgi)

USE types
IMPLICIT NONE

INTEGER                                         :: cur_index
INTEGER                                         :: pack_index, cur_pgi
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: pos, corner, width

DOUBLE PRECISION                                :: diff1, diff2

pos = package(pack_index)%pos
corner = dyn_cell(cur_pgi)%corner
width = dyn_cell(cur_pgi)%width


 diff1 = abs(pos(cur_index) - corner(cur_index))
 diff2 = abs(pos(cur_index) - corner(cur_index) - width(cur_index))
write(*,*) 'correction_propagation: diff1 = ', diff1, ' diff2 = ', diff2

write(*,*) 'correction_propagation: corrected index: ', cur_index

IF(diff1 < diff2) THEN
 package(pack_index)%pos(cur_index) = corner(cur_index)
 write(*,*) 'correction_propagation: new pos = corner(cur_index)'
ELSE IF(diff2 >= diff1) THEN
 package(pack_index)%pos(cur_index) = corner(cur_index) + width(cur_index)
 write(*,*) 'correction_propagation: new pos = corner(cur_index) + width(cur_index)'
END IF

















END SUBROUTINE correction_propagation
