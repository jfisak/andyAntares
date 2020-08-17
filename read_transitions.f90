SUBROUTINE read_transitions(el_index, lowerion, upperion, transition_type, filename)


 USE types

 IMPLICIT NONE    

 ! input data
 INTEGER                        :: element, lowerion, upperion, transition_type, el_index
 INTEGER                        :: n_ions
 CHARACTER (LEN=30)             :: filename
 CHARACTER (LEN=300)            :: line
 INTEGER                        :: jint
 DOUBLE PRECISION               :: jdble
 REAL                           :: jfloat
 ! used constants
 DOUBLE PRECISION               :: oconstant
 ! loop variables
 INTEGER                        :: I, J
 INTEGER, ALLOCATABLE           :: ntrans(:)
 TYPE(line_list), ALLOCATABLE   :: pom(:)
 INTEGER                        :: new_size, old_size
 INTEGER                        :: current_index, current_ntrans
 INTEGER                        :: trans
 INTEGER                        :: ion_index
 LOGICAL                        :: found_low_conf, found_up_conf
 ! reading from the file
 INTEGER                        :: at_number
 INTEGER                        :: reading_transitions, n_levels, tot_ntrans
 INTEGER                        :: current_element, current_ion
 INTEGER                        :: kindex, low_level, up_level
 ! actual lower and upper index
 INTEGER                        :: act_lower, act_upper, low_vsplit, up_vsplit
 ! counter of not included transitions
 INTEGER                        :: n_not_included
 CHARACTER (LEN=16)             :: low_conf
 CHARACTER (LEN=16)             :: up_conf
 DOUBLE PRECISION               :: A, col_str, l_freq
 INTEGER                        :: electron_number
 DOUBLE PRECISION               :: ee_lc, ee_uc
 INTEGER                        :: act_lc, act_uc
 DOUBLE PRECISION               :: g_lower
 DOUBLE PRECISION               :: deltaE

992 format(I7, I3, I3, I5, I5, I4, I4, Tr1, A16, A16, e9.2, e10.2, e10.2, f5.1, f5.1)

orbitals_nl = .FALSE.

! calculate the constant for the oscilator strength calculation
oconstant = (me_g * light_speed ** 3)/(8.D0 * pi ** 2 * e_charge**2)
n_ions = 0
! element: atomic number, could be different from el_index
element = elements(el_index)%atom_number
simpleTrans = .FALSE.
OPEN (UNIT=9, status='old', FILE=filename)
IF(transition_type == 9) THEN
 simpleTrans = .TRUE.
 transition_type = 2
END IF


SELECT CASE(transition_type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! standard atomic data (26.03.2016)
! **0** we ignore rows starting *
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)
 STOP 'reading_transitions: the standard type of data is no longer supported'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data from the Opacity project
! **2**
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 ! computing number of transitions for the given ions
 ! initialization of number of transitions
 orbitals_nl = .TRUE.
 ALLOCATE(ntrans(upperion - lowerion + 1))
 DO I = lowerion, upperion
  ion_index = I - lowerion + 1
  ntrans(ion_index) = 0
 END DO
 tot_ntrans = 0
 DO
  READ(9,992, IOSTAT = reading_transitions) jint, at_number, electron_number,&
   jint, low_vsplit, up_vsplit, jint, low_conf, up_conf, jdble, jdble
  ! write(*,*) 'read_transitions: low_conf = ', TRIM(low_conf), ' up_conf = ', TRIM(up_conf)
  ! write(*,*) 'read_transitions: at_number = ', at_number
  ! READ(9,'(A)',IOSTAT = reading_transitions) line
  !print*, line
  IF(reading_transitions /= 0) EXIT
  ! IF(line(1:1) == '*') CYCLE
  CALL find_element_index(at_number,el_index)
  current_ion = element - electron_number + 1
  ! write(99,*) 'read_transitions: current_ion = ', current_ion
  found_low_conf = .FALSE.
  found_up_conf = .FALSE.
  n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)
  DO J = 1, n_levels
   ! write(*,*) 'read_transitions: low_conf = ', low_conf, ' up_conf = ', up_conf, &
   !  ' elconf = ', elements(el_index)%ions(current_ion)%levels(J)%elconf
   IF(low_conf .EQ. elements(el_index)%ions(current_ion)%levels(J)%elconf &
     .AND. low_vsplit == elements(el_index)%ions(current_ion)%levels(J)%vsplit) THEN
    found_low_conf = .TRUE.
    ! write(*,*) 'found low_conf ', low_conf
    ! write(*,*) 'found electron configuration...'
   END IF
   IF(TRIM(up_conf) .EQ. TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf) &
     .AND. up_vsplit == elements(el_index)%ions(current_ion)%levels(J)%vsplit) THEN
    found_up_conf = .TRUE.
    ! write(*,*) 'found up_conf ', up_conf
    ! write(*,*) 'found electron configuration...'
   END IF
   !IF(found_up_conf .EQV. .TRUE. .AND. found_low_conf .EQV. .TRUE.) EXIT
   IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) EXIT
  END DO
  !write(99,*) 'low_conf = ', low_conf, ' up_conf = ', up_conf, found_up_conf, found_low_conf
  IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) THEN
   ! write(*,*) 'reading_transitions: el = ', el_index, ' ion = ', current_ion,&
   !  ' line ', low_conf, ' -> ', up_conf, ' was accepted'
   ion_index = current_ion - lowerion + 1
   ntrans(current_ion) = ntrans(current_ion) + 1
   tot_ntrans = tot_ntrans + 1
  END IF
  ! write(*,*) 'tot_ntrans = ', tot_ntrans
 END DO
 ! we have to treat the linelist array
 ! 1.) if it is allocated we add to this array tot_ntrans new fields
 ! 2.) we allocate a new array with dimension tot_ntrans otherwise
 ! ntransitions is a "pointer" (and a global variable) which points to the last saved transition
 IF(ALLOCATED(linelist)) THEN
  old_size = SIZE(linelist)
  new_size = old_size + tot_ntrans
  ALLOCATE(pom(old_size))
  pom(1:old_size) = linelist(1:old_size)
  DEALLOCATE(linelist)
  ALLOCATE(linelist(new_size))
  linelist(1:old_size) = pom(1:old_size)
  DEALLOCATE(pom)
 ELSE
  ALLOCATE(linelist(tot_ntrans))
  ntransitions = 0
 END IF
 REWIND(9)
 ! we have \sum_i ntrans_i additional transitions into the linelist array
 ! now we will read the given data
 DO I = lowerion, upperion
  current_index = I - lowerion + 1
  current_ntrans = ntrans(current_index)
  ! trans is a number of read transitions 
  ! because we have tu assume, that we can comment lines in the middle of the file
  trans = 0
  n_not_included = 0
  DO
   READ(9,992, IOSTAT = reading_transitions) kindex, current_element, electron_number,&
   low_vsplit, up_vsplit, low_level, up_level, low_conf, up_conf, col_str, A, l_freq
   ! write(*,*) 'read_transitions: col_str = ', col_str, ' A = ', A, &
   !  ' l_freq = ', l_freq
   ! READ(9,'(A)',IOSTAT = reading_transitions) line
   ! write(99,*) line
   IF(reading_transitions /= 0) EXIT
   ! IF(line(1:1) == '*') CYCLE
   ! READ(line,*) kindex, current_element, electron_number, junk, junk, low_level, up_level, &
   !      low_conf, up_conf, col_str, A, l_freq
   current_ion = element - electron_number + 1
   ! write(99,*) 'I = ', I, 'i_conf = ', low_conf, 'j_conf = ', up_conf, &
   ! 'l_index = ', elements(el_index)%ions(current_ion)%levels(I)%l_index

   ! firstly we will check out if this is a transition which levels were really included
   ! we do not want to save lines with no probability
   IF(A == 0.D0) THEN
    !write(99,*) 'A is equal to zero, cycling...'
    CYCLE
   END IF
   ! write(*,*) 'I = ', I, 'i_conf = ', low_conf, 'j_conf = ', up_conf
   ! because we could find all included transitions before the 
   IF(ntransitions == SIZE(linelist)) EXIT
   found_low_conf = .FALSE.
   found_up_conf = .FALSE.
   n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)
   DO J = 1, n_levels
    IF(TRIM(low_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf) &
     .AND. low_vsplit == elements(el_index)%ions(current_ion)%levels(J)%vsplit) THEN
     found_low_conf = .TRUE.
     act_lc = J
     !write(99,*) 'found electron configuration...'
     ! IF(col_str >= 0) THEN
     !  act_upper = elements(el_index)%ions(current_ion)%levels(J)%l_index
     ! ELSE
     !  act_lower = elements(el_index)%ions(current_ion)%levels(J)%l_index
     ! END IF 
     !write(99,*) 'l_index = ', elements(el_index)%ions(current_ion)%levels(J)%l_index
    END IF
    IF(TRIM(up_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf) &
     .AND. up_vsplit == elements(el_index)%ions(current_ion)%levels(J)%vsplit) THEN
     found_up_conf = .TRUE.
     act_uc = J
     !write(99,*) 'found electron configuration...'
     ! IF(col_str >= 0) THEN
     !  act_lower = elements(el_index)%ions(current_ion)%levels(J)%l_index
     ! ELSE
     !  act_upper = elements(el_index)%ions(current_ion)%levels(J)%l_index
     ! END IF
     !write(99,*) 'l_index = ', elements(el_index)%ions(current_ion)%levels(J)%l_index
    END IF
    IF((found_low_conf .EQV. .TRUE.) .AND. (found_up_conf .EQV. .TRUE.)) THEN
     ee_lc = elements(el_index)%ions(current_ion)%levels(act_lc)%exci_energy
     ee_uc = elements(el_index)%ions(current_ion)%levels(act_uc)%exci_energy
     IF(ee_lc < ee_uc) THEN
      act_lower = act_lc
      act_upper = act_uc
     ELSE IF(ee_uc < ee_lc) THEN
      act_lower = act_uc
      act_upper = act_lc
     ELSE
      write(*,*) 'read_transitions: act_lc = ', act_lc, ' act_uc = ', act_uc
      write(*,*) 'conf_l = ', TRIM(low_conf), ' conf_u = ', up_conf
      write(*,*) 'ee_lc = ', ee_lc, ' ee_uc = ', ee_uc
      STOP 'ee_lc == ee_uc'
     END IF
    END IF
   END DO
   IF((found_low_conf .EQV. .TRUE.) .AND. (found_up_conf .EQV. .TRUE.)) THEN
    ! write(*,*) 'read_transitions: element: ', element, ' ion = ', ion_index, &
    !  ' line from ', low_conf, ' to ', up_conf, 'lc = ', act_lc, 'uc = ', act_uc, ' was included...'
   END IF
   IF((found_low_conf .EQV. .FALSE.) .OR. (found_up_conf .EQV. .FALSE.)) THEN
    ! this configuration will not be taken into account and we will read the next line
    write(*,*) 'element: ', element, ' ion = ', ion_index, ' line from ', low_conf, ' to ', up_conf, &
      'low_vsplit = ', low_vsplit, ' up_vsplit = ', up_vsplit, '  was not included...'
    n_not_included = n_not_included + 1
    CYCLE
   ELSE
    ntransitions = ntransitions + 1
    ! write(*,*) 'reading_transitions: ntransitions = ', ntransitions
    linelist(ntransitions)%lower = act_lower
    linelist(ntransitions)%upper = act_upper
   END IF
   IF((current_element /= element) .OR. (current_ion < lowerion) &
     .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC TRANSITIONS...'
   !write(99,*) 'indexe = ', el_index, ' indexi = ', current_ion, ' ntransitions = ', ntransitions
   g_lower = elements(el_index)%ions(current_ion)%levels(act_lower)%stat_waight
   linelist(ntransitions)%indexe = el_index
   linelist(ntransitions)%indexi = current_ion
   ! linelist(ntransitions)%freq = 1.E+8 * light_speed / l_freq
   deltaE = elements(el_index)%ions(current_ion)%levels(act_upper)%exci_energy - &
    elements(el_index)%ions(current_ion)%levels(act_lower)%exci_energy
   linelist(ntransitions)%freq =  deltaE / h
   ! write(*,*) 'read_transitions: lambda = ', light_speed / linelist(ntransitions)%freq * 1.E8
   IF(simpleTrans) THEN
    linelist(ntransitions)%A_ul = abs(A)
    linelist(ntransitions)%f_ul = abs(col_str)
   ELSE
    linelist(ntransitions)%A_ul = abs(A) / g_lower
    linelist(ntransitions)%f_ul = abs(col_str) / g_lower
   END IF
   linelist(ntransitions)%n_int = 0
   ! write(99,*) 'line: ', ntransitions, ' el = ', el_index, ' ion = ', current_ion,&
   !  ' lower level energy = ', &
   !  elements(el_index)%ions(ion_index)%levels(linelist(ntransitions)%lower)%exci_energy / e_v,& 
   !  elements(el_index)%ions(ion_index)%levels(linelist(ntransitions)%upper)%exci_energy / e_v
   ! we have a transition between two atomic levels and we have to connect
   ! this transition with already read levels from another file, we do this
   ! connection via electron configuration in the form of string
   trans = trans + 1
   ! write(99,*) 'read_transitions: trans = ', trans, ' ntrans(I) = ', ntrans(I)
   IF(trans == ntrans(I)) EXIT
  END DO
  ! write(99,*) 'read_transitions: el = ', el_index, ' ion = ', current_ion, n_not_included, 'lines were not included'
 END DO
 !STOP 'read_transitions: testing'
 ! if we did now use every transition in the file we will reallocate the array
 ! linelist so it will not be so large
 IF(ntransitions < SIZE(linelist)) THEN
  ALLOCATE(pom(SIZE(linelist)))
  pom(1:ntransitions) = linelist(1:ntransitions)
  DEALLOCATE(linelist)
  ALLOCATE(linelist(ntransitions))
  linelist(1:ntransitions) = pom(1:ntransitions)
  DEALLOCATE(pom)
 END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! data from the Tardis code
! **8**
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(8)
 ! computing number of transitions for the given ions
 ! initialization of number of transitions
 write(*,*) 'sbr read_transitions...'
 ALLOCATE(ntrans(upperion - lowerion + 1))
 DO I = lowerion, upperion
  ion_index = I
  ntrans(ion_index) = 0
 END DO
 tot_ntrans = 0
 ! calculation a number of transition which will be included
 DO ! main loop
  READ(9,'(A)',IOSTAT = reading_transitions) line
  ! write(*,*) 'read_transitions: line = ', line
  IF(reading_transitions /= 0) EXIT
  IF(line(1:1) == '*') CYCLE
  READ(line,*) at_number, electron_number,&
   low_conf,  up_conf, jint, jfloat, jfloat, jfloat,  jdble, jdble, jdble, jdble, jfloat
  ! write(*,*) 'read_transitions: low_conf = ', TRIM(low_conf), ' up_conf = ', TRIM(up_conf)
  write(*,*) 'read_transitions: electron_number = ', electron_number
  !print*, line
  write(*,*) 'read_transitions: reading_transitions = ', reading_transitions
  CALL find_element_index(at_number,el_index)
  current_ion = electron_number + 1
  ! write(99,*) 'read_transitions: current_ion = ', current_ion
  found_low_conf = .FALSE.
  found_up_conf = .FALSE.
  n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)
  DO J = 1, n_levels
   ! write(*,*) 'read_transitions: low_conf = ', low_conf, ' up_conf = ', up_conf, &
   !  ' elconf = ', elements(el_index)%ions(current_ion)%levels(J)%elconf
   IF(TRIM(low_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf)) THEN
    found_low_conf = .TRUE.
    ! write(*,*) 'found low_conf ', low_conf
    ! write(*,*) 'found electron configuration...'
   END IF
   IF(TRIM(up_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf)) THEN
    found_up_conf = .TRUE.
    ! write(*,*) 'found up_conf ', up_conf
    ! write(*,*) 'found electron configuration...'
   END IF
   !IF(found_up_conf .EQV. .TRUE. .AND. found_low_conf .EQV. .TRUE.) EXIT
   IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) EXIT
  END DO
  !write(99,*) 'low_conf = ', low_conf, ' up_conf = ', up_conf, found_up_conf, found_low_conf
  IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) THEN
   ! write(*,*) 'reading_transitions: el = ', el_index, ' ion = ', current_ion,&
   !  ' line ', low_conf, ' -> ', up_conf, ' was accepted'
   ion_index = current_ion - lowerion + 1
   ntrans(current_ion) = ntrans(current_ion) + 1
   tot_ntrans = tot_ntrans + 1
  END IF
  ! write(*,*) 'tot_ntrans = ', tot_ntrans
 END DO ! end main loop

 ! we have to treat the linelist array
 ! 1.) if it is allocated we add to this array tot_ntrans new fields
 ! 2.) we allocate a new array with dimension tot_ntrans otherwise
 ! ntransitions is a "pointer" (and a global variable) which points to the last saved transition
 IF(ALLOCATED(linelist)) THEN
  old_size = SIZE(linelist)
  new_size = old_size + tot_ntrans
  ALLOCATE(pom(old_size))
  pom(1:old_size) = linelist(1:old_size)
  DEALLOCATE(linelist)
  ALLOCATE(linelist(new_size))
  linelist(1:old_size) = pom(1:old_size)
  DEALLOCATE(pom)
 ELSE
  ALLOCATE(linelist(tot_ntrans))
  ntransitions = 0
 END IF
 REWIND(9)

 ! we have \sum_i ntrans_i additional transitions into the linelist array
 ! now we will read the given data
 ! DO I = lowerion, upperion ! do 01 over ions
 !  current_index = I - lowerion + 1
 !  current_ntrans = ntrans(current_index)
  ! trans is a number of read transitions 
  ! because we have tu assume, that we can comment lines in the middle of the file
  trans = 0
  n_not_included = 0
  DO ! do 02 reading line by line
   READ(9,'(A)',IOSTAT = reading_transitions) line
   write(*,*) 'read_transitions: line = ', line
   IF(reading_transitions /= 0) EXIT
   write(*,*) 'read_transitions: reading_transitions = ', reading_transitions
   IF(line(1:1) == '*') CYCLE
   READ(line, *) current_element, electron_number,&
    low_conf,  up_conf, jint, jdble, col_str, jdble,  jdble, jdble, jdble, A, jdble
   ! write(*,*) 'read_transitions: col_str = ', col_str, ' A = ', A, &
   !  ' l_freq = ', l_freq
   ! READ(9,'(A)',IOSTAT = reading_transitions) line
   ! write(99,*) line
   IF(reading_transitions /= 0) EXIT
   ! IF(line(1:1) == '*') CYCLE
   ! READ(line,*) kindex, current_element, electron_number, junk, junk, low_level, up_level, &
   !      low_conf, up_conf, col_str, A, l_freq
   current_ion = electron_number + 1
   ! write(99,*) 'I = ', I, 'i_conf = ', low_conf, 'j_conf = ', up_conf, &
   ! 'l_index = ', elements(el_index)%ions(current_ion)%levels(I)%l_index

   ! firstly we will check out if this is a transition which levels were really included
   ! we do not want to save lines with no probability
   IF(A == 0.D0) THEN
    !write(99,*) 'A is equal to zero, cycling...'
    CYCLE
   END IF
   ! write(*,*) 'I = ', I, 'i_conf = ', low_conf, 'j_conf = ', up_conf
   ! because we could find all included transitions before the 
   IF(ntransitions == SIZE(linelist)) EXIT
   found_low_conf = .FALSE.
   found_up_conf = .FALSE.
   n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)
   DO J = 1, n_levels ! loop 03 over levels
    write(*,*) 'read_transitions: low_conf = ', low_conf, ' up_conf = ', up_conf, &
     ' elconf = ', elements(el_index)%ions(current_ion)%levels(J)%elconf
    IF(TRIM(low_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf)) THEN
     found_low_conf = .TRUE.
     act_lc = J
     ! write(*,*) 'found low_conf ', low_conf
     ! write(*,*) 'found electron configuration...'
    END IF
    IF(TRIM(up_conf) == TRIM(elements(el_index)%ions(current_ion)%levels(J)%elconf)) THEN
     found_up_conf = .TRUE.
     act_uc = J
     ! write(*,*) 'found up_conf ', up_conf
     ! write(*,*) 'found electron configuration...'
    END IF
    !IF(found_up_conf .EQV. .TRUE. .AND. found_low_conf .EQV. .TRUE.) EXIT
    IF((found_up_conf .EQV. .TRUE.) .AND. (found_low_conf .EQV. .TRUE.)) EXIT
    IF((found_low_conf .EQV. .TRUE.) .AND. (found_up_conf .EQV. .TRUE.)) THEN
     ee_lc = elements(el_index)%ions(current_ion)%levels(act_lc)%exci_energy
     ee_uc = elements(el_index)%ions(current_ion)%levels(act_uc)%exci_energy
     IF(ee_lc < ee_uc) THEN
      act_lower = act_lc
      act_upper = act_uc
     ELSE IF(ee_uc < ee_lc) THEN
      act_lower = act_uc
      act_upper = act_lc
     ELSE
      write(*,*) 'read_transitions: act_lc = ', act_lc, ' act_uc = ', act_uc
      write(*,*) 'conf_l = ', TRIM(low_conf), ' conf_u = ', up_conf
      write(*,*) 'ee_lc = ', ee_lc, ' ee_uc = ', ee_uc
      STOP 'ee_lc == ee_uc'
     END IF
    END IF
   END DO ! end loop 03 over levels
   IF((found_low_conf .EQV. .TRUE.) .AND. (found_up_conf .EQV. .TRUE.)) THEN
     write(*,*) 'read_transitions: element: ', element, ' ion = ', ion_index, &
      ' line from ', low_conf, ' to ', up_conf, 'lc = ', act_lc, 'uc = ', act_uc, ' was included...'
    ee_lc = elements(el_index)%ions(current_ion)%levels(act_lc)%exci_energy
    ee_uc = elements(el_index)%ions(current_ion)%levels(act_uc)%exci_energy
    IF(ee_lc < ee_uc) THEN
     act_lower = act_lc
     act_upper = act_uc
    ELSE IF(ee_uc < ee_lc) THEN
     act_lower = act_uc
     act_upper = act_lc
    ELSE
     STOP 'read_transitions: ee_lc == ee_uc'
    END IF
   END IF
   IF((found_low_conf .EQV. .FALSE.) .OR. (found_up_conf .EQV. .FALSE.)) THEN
    ! this configuration will not be taken into account and we will read the next line
    write(*,*) 'element: ', element, ' ion = ', ion_index, ' line from ', low_conf, ' to ', up_conf, &
      'low_vsplit = ', low_vsplit, ' up_vsplit = ', up_vsplit, '  was not included...'
    n_not_included = n_not_included + 1
    CYCLE
   ELSE
    ntransitions = ntransitions + 1
    ! write(*,*) 'reading_transitions: ntransitions = ', ntransitions
    linelist(ntransitions)%lower = act_lower
    linelist(ntransitions)%upper = act_upper
   END IF
   IF((current_element /= element) .OR. (current_ion < lowerion) &
     .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC TRANSITIONS...'
   !write(99,*) 'indexe = ', el_index, ' indexi = ', current_ion, ' ntransitions = ', ntransitions
   g_lower = elements(el_index)%ions(current_ion)%levels(act_lower)%stat_waight
   linelist(ntransitions)%indexe = el_index
   linelist(ntransitions)%indexi = current_ion
   ! linelist(ntransitions)%freq = 1.E+8 * light_speed / l_freq
   deltaE = elements(el_index)%ions(current_ion)%levels(act_upper)%exci_energy - &
    elements(el_index)%ions(current_ion)%levels(act_lower)%exci_energy
   linelist(ntransitions)%freq =  deltaE / h
   linelist(ntransitions)%A_ul = A
   linelist(ntransitions)%f_ul = col_str
   ! write(*,*) 'read_transitions: lambda = ', light_speed / linelist(ntransitions)%freq * 1.E8
   linelist(ntransitions)%n_int = 0
   ! write(99,*) 'line: ', ntransitions, ' el = ', el_index, ' ion = ', current_ion,&
   !  ' lower level energy = ', &
   !  elements(el_index)%ions(ion_index)%levels(linelist(ntransitions)%lower)%exci_energy / e_v,& 
   !  elements(el_index)%ions(ion_index)%levels(linelist(ntransitions)%upper)%exci_energy / e_v

   ! we have a transition between two atomic levels and we have to connect
   ! this transition with already read levels from another file, we do this
   ! connection via electron configuration in the form of string
   trans = trans + 1
   ! write(99,*) 'read_transitions: trans = ', trans, ' ntrans(I) = ', ntrans(I)
  END DO
  ! write(99,*) 'read_transitions: el = ', el_index, ' ion = ', current_ion, n_not_included, 'lines were not included'
 ! END DO
 !STOP 'read_transitions: testing'
 ! if we did now use every transition in the file we will reallocate the array
 ! linelist so it will not be so large
 IF(ntransitions < SIZE(linelist)) THEN
  ALLOCATE(pom(SIZE(linelist)))
  pom(1:ntransitions) = linelist(1:ntransitions)
  DEALLOCATE(linelist)
  ALLOCATE(linelist(ntransitions))
  linelist(1:ntransitions) = pom(1:ntransitions)
  DEALLOCATE(pom)
 END IF
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! **DEFAULT**
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE DEFAULT
  write(*,*) 'no transitions found'
  STOP
 END SELECT
CLOSE(9) 

END SUBROUTINE read_transitions
