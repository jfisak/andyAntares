SUBROUTINE read_singleline()

USE types
IMPLICIT NONE

INTEGER                         :: kindex, current_element, electron_number
INTEGER                         :: current_ion
INTEGER                         :: ind_J
INTEGER                         :: low_vsplit, up_vsplit, low_level, up_level
CHARACTER (LEN=16)              :: low_conf
CHARACTER (LEN=16)              :: up_conf
DOUBLE PRECISION                :: col_str, trans_A, l_freq
LOGICAL                         :: found_low_conf, found_up_conf

INTEGER                         :: el_index
INTEGER                         :: ion_index, g_lower
INTEGER                         :: n_levels, tot_ntrans, reading_transitions


992 format(I7, I3, I3, I5, I5, I4, I4, Tr1, A16, A16, e9.2, e10.2, e10.2, f5.1, f5.1)

tot_ntrans = 1
ALLOCATE(linelist(tot_ntrans))

OPEN(9)
 READ(9,992, IOSTAT = reading_transitions) kindex, current_element, electron_number,&
  low_vsplit, up_vsplit, low_level, up_level, low_conf, up_conf, col_str, trans_A, l_freq
CLOSE(9)

CALL find_element_index(current_element, el_index)
current_ion = current_element - electron_number + 1

found_low_conf = .FALSE.
found_up_conf = .FALSE.

n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)
DO ind_J = 1, n_levels
 IF(TRIM(low_conf) == elements(el_index)%ions(current_ion)%levels(ind_J)%elconf &
   .AND. low_vsplit == elements(el_index)%ions(current_ion)%levels(ind_J)%vsplit) THEN
  found_low_conf = .TRUE.
 END IF
 IF(TRIM(up_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(ind_J)%elconf) &
   .AND. up_vsplit == elements(el_index)%ions(current_ion)%levels(ind_J)%vsplit) THEN
  found_up_conf = .TRUE.
 END IF
 IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) EXIT
END DO
IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) THEN
 ion_index = current_ion
END IF

IF(simpleTrans) THEN
 linelist(ntransitions)%A_ul = abs(trans_A)
 linelist(ntransitions)%f_lu = abs(col_str)
ELSE
 linelist(ntransitions)%A_ul = abs(trans_A) / g_lower
 linelist(ntransitions)%f_lu = abs(col_str) / g_lower
END IF







END SUBROUTINE read_singleline
