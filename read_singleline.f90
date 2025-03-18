SUBROUTINE read_singleline()

USE types
USE constants
IMPLICIT NONE

INTEGER                         :: kindex, current_element, electron_number
INTEGER                         :: current_ion
INTEGER                         :: ind_J
INTEGER                         :: low_vsplit, up_vsplit, low_level, up_level
CHARACTER (LEN=16)              :: low_conf
CHARACTER (LEN=16)              :: up_conf
DOUBLE PRECISION                :: col_str, trans_A, l_freq, deltaE
LOGICAL                         :: found_low_conf, found_up_conf

INTEGER                         :: el_index
INTEGER                         :: n_levels, tot_ntrans, reading_transitions
INTEGER                         :: act_lc, act_uc, act_lower, act_upper
DOUBLE PRECISION                :: ee_lc, ee_uc

992 format(I7, I3, I3, I5, I5, I4, I4, Tr1, A16, A16, e9.2, e10.2, e10.2, f5.1, f5.1)

tot_ntrans = 1
ntransitions = 1
ALLOCATE(linelist(tot_ntrans))

OPEN(UNIT=9, STATUS='old', FILE=singleline_file)
 READ(9,992, IOSTAT = reading_transitions) kindex, current_element, electron_number,&
  low_vsplit, up_vsplit, low_level, up_level, low_conf, up_conf, col_str, trans_A, l_freq
  write(*,*) 'read_singleline: kindex = ', kindex, ' current_element = ', current_element
CLOSE(9)

write(*,*) 'read_singleline: current_element = ', current_element
CALL find_element_index(current_element, el_index)
current_ion = current_element - electron_number + 1

found_low_conf = .FALSE.
found_up_conf = .FALSE.

n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)

DO ind_J = 1, n_levels ! loop 03 over levels
 ! write(*,*) 'read_transitions: low_conf = ', low_conf, ' up_conf = ', up_conf, &
 !  ' elconf = ', elements(el_index)%ions(current_ion)%levels(ind_J)%elconf
 IF(TRIM(low_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(ind_J)%elconf)) THEN
  found_low_conf = .TRUE.
  act_lc = ind_J
  ! write(*,*) 'found low_conf ', low_conf
  ! write(*,*) 'found electron configuration...'
 END IF
 IF(TRIM(up_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(ind_J)%elconf)) THEN
  found_up_conf = .TRUE.
  act_uc = ind_J
  ! write(*,*) 'found up_conf ', up_conf
  ! write(*,*) 'found electron configuration...'
 END IF
 IF((found_low_conf .EQV. .TRUE.) .AND. (found_up_conf .EQV. .TRUE.)) THEN
  ee_lc = elements(el_index)%ions(current_ion)%levels(act_lc)%exci_energy
  ee_uc = elements(el_index)%ions(current_ion)%levels(act_uc)%exci_energy
  IF(ee_lc < ee_uc) THEN
   act_lower = act_lc
   act_upper = act_uc
   EXIT
  ELSE IF(ee_uc < ee_lc) THEN
   act_lower = act_uc
   act_upper = act_lc
   EXIT
  ELSE
   ! write(*,*) 'read_transitions: act_lc = ', act_lc, ' act_uc = ', act_uc
   ! write(*,*) 'conf_l = ', TRIM(low_conf), ' conf_u = ', up_conf
   ! write(*,*) 'ee_lc = ', ee_lc, ' ee_uc = ', ee_uc
   STOP 'ee_lc == ee_uc'
  END IF
 END IF
END DO ! end loop 03 over levels

linelist(tot_ntrans)%lower = act_lower
linelist(tot_ntrans)%upper = act_upper

linelist(tot_ntrans)%indexe = el_index
linelist(tot_ntrans)%indexi = current_ion
deltaE = elements(el_index)%ions(current_ion)%levels(act_upper)%exci_energy - &
 elements(el_index)%ions(current_ion)%levels(act_lower)%exci_energy
linelist(tot_ntrans)%freq =  deltaE / const_h
linelist(tot_ntrans)%A_ul = trans_A
linelist(tot_ntrans)%f_lu = col_str
linelist(tot_ntrans)%n_int = 0






END SUBROUTINE read_singleline
