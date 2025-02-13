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
INTEGER                                         :: pack_index, cur_pgi, new_pgi
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: corner, width, cur_pos

DOUBLE PRECISION                                :: diff1, diff2

cur_pos = package(pack_index)%pos
corner = dyn_cell(cur_pgi)%corner
width = dyn_cell(cur_pgi)%width

IF(debug == 2) THEN
 write(*,*) 'correction_propagation: cell starting = ', corner
 write(*,*) 'correction_propagation: packet pos = ', cur_pos
 write(*,*) 'correction_propagation: cell ending = ', (corner + width)
END IF

diff1 = abs(cur_pos(cur_index) - corner(cur_index))
diff2 = abs(cur_pos(cur_index) - corner(cur_index) - width(cur_index))
IF(debug == 2) THEN
 write(*,*) 'correction_propagation: diff1 = ', diff1, ' diff2 = ', diff2
 
 write(*,*) 'correction_propagation: corrected index: ', cur_index
END IF

IF(diff1 < diff2) THEN
 package(pack_index)%pos(cur_index) = corner(cur_index)
 IF(debug == 2) THEN
  write(*,*) 'correction_propagation: new pos = corner(cur_index)'
 END IF
ELSE IF(diff1 >= diff2) THEN
 package(pack_index)%pos(cur_index) = corner(cur_index) + width(cur_index)
 IF(debug == 2) THEN
  write(*,*) 'correction_propagation: new pos = corner(cur_index) + width(cur_index)'
 END IF
END IF


CALL find_dyn_cell1(cur_pos, new_pgi)
package(pack_index)%cell_numb = new_pgi
IF(debug == 2) THEN
 write(*,*) 'correction_propagation: cell starting = ', corner
 write(*,*) 'correction_propagation: packet pos = ', package(pack_index)%pos
 write(*,*) 'correction_propagation: cell ending = ', (corner + width)
END IF

END SUBROUTINE correction_propagation
