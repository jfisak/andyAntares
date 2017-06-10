SUBROUTINE read_transitions(el_index, lowerion, upperion, transition_type, filename)


 USE types

 IMPLICIT NONE    

 ! input data
 INTEGER                        :: element, lowerion, upperion, transition_type, el_index
 INTEGER                        :: n_ions
 CHARACTER (LEN=20)             :: filename, junk
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
 CHARACTER (LEN=200)            :: line
 INTEGER                        :: ios, reading_transitions, linereading, n_levels, tot_ntrans
 INTEGER                        :: n_transitions, curr_n_tran, n_line
 INTEGER                        :: current_element, current_ion
 INTEGER                        :: kindex, low_level, up_level
 CHARACTER (LEN=15)              :: low_conf, up_conf
 DOUBLE PRECISION               :: A, col_str, l_freq


! calculate the constant for the oscilator strength calculation
oconstant = (me_g * light_speed ** 3)/(8.D0 * pi ** 2 * e_charge**2)
n_ions = 0
element = elements(el_index)%atom_number
OPEN (UNIT=9, status='old', FILE=filename)
SELECT CASE(transition_type)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! standard atomic data (26.03.2016)
! **0** we ignore rows starting *
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(0)
 ! only if the ray linelist is allocated we can read important physical quantities from file
 IF(ALLOCATED(linelist)) THEN
 DO
  READ(9,'(A)',IOSTAT=ios) line
  IF (ios /= 0) EXIT
  IF ( INDEX(line, '*') /=0) CYCLE
  READ(line,*) current_element, current_ion, n_transitions
  IF(n_transitions == 0) CYCLE
  ! do we have the right file?
  IF((current_element /= element) .OR. (current_ion < lowerion) &
    .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC TRANSITIONS...'
  n_ions = elements(el_index)%nions
  !print*, 'number of ions n_ions = ', n_ions
   DO I = 1, n_transitions
    READ(9,'(A)',IOSTAT=linereading) line
     IF (linereading /= 0) THEN
     print*, 'the subroutine read_transitions:'
     STOP 'ERROR: NO VALID TRANSITIONS...'
    END IF
    READ(line,*) kindex, low_level, up_level, A, col_str
    !n_line = linelist(ntransitions)%indexe + 1
    !linelist(ntransitions)%indexe = n_line
    n_line = ntransitions + 1
    ntransitions = n_line
    ! save data to the linelist
    linelist(n_line)%indexe = current_element
    linelist(n_line)%indexi = current_ion
    linelist(n_line)%lower = low_level
    linelist(n_line)%upper = up_level
    linelist(n_line)%freq = (elements(el_index)%ions(current_ion)%levels(up_level)%exci_energy -&
                               elements(el_index)%ions(current_ion)%levels(low_level)%exci_energy) / h
    linelist(n_line)%A_ul = A
    linelist(n_line)%f_ul = oconstant * (elements(el_index)%ions(current_ion)%levels(up_level)%stat_waight / &
                       elements(el_index)%ions(current_ion)%levels(low_level)%stat_waight) *                 &
                       (linelist(n_line)%A_ul / linelist(n_line)%freq ** 2)
    !print*, n_line, ': ', 'indexe: ', linelist(n_line)%indexe, linelist(n_line)%indexi, linelist(n_line)%lower, &
    !                   linelist(n_line)%upper, 'f = ', linelist(n_line)%freq, linelist(n_line)%A_ul, linelist(n_line)%f_ul
      write(20,*) linelist(n_line)%freq, linelist(n_line)%f_ul
   END DO
   n_ions = n_ions + 1
   IF(n_ions == (upperion - lowerion + 1)) EXIT
 END DO
 ! if we do not have the ray linelist allocated, we have
 ! to calculate number of possible transitions
 ELSE
 DO   
  READ(9,'(A)',IOSTAT=ios) line
  !print*, line
  IF (ios /= 0) EXIT
  IF ( INDEX(line, '*') /= 0) CYCLE
  READ(line,*) junk, junk, n_transitions
  ntransitions = ntransitions + n_transitions
  IF(n_transitions == 0) CYCLE
  ! loop which will do nothing
  DO I=1,n_transitions
   READ(9,*)
  END DO
 END DO
 END IF
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
 ALLOCATE(ntrans(upperion - lowerion + 1))
 DO I = lowerion, upperion
  ion_index = I - lowerion + 1
  ntrans(ion_index) = 0
 END DO
 tot_ntrans = 0
 DO
  READ(9,'(A)',IOSTAT = reading_transitions) line
  !print*, line
  IF(reading_transitions /= 0) EXIT
  IF(line(1:1) == '*') CYCLE
  READ(line,*) junk, junk, current_ion, junk, junk, junk, junk, junk, junk, junk, junk, junk 
  ion_index = current_ion - lowerion + 1
  ntrans(current_ion) = ntrans(current_ion) + 1
  tot_ntrans = tot_ntrans + 1
  !print*, 'tot_ntrans = ', tot_ntrans
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
  DO
   READ(9,'(A)',IOSTAT = reading_transitions) line
   !print*, line
   IF(reading_transitions /= 0) EXIT
   IF(line(1:1) == '*') CYCLE
   READ(line,*) kindex, current_element, current_ion, junk, junk, low_level, up_level, &
        low_conf, up_conf, col_str, A, l_freq
   !print*, 'I = ', I, 'i_conf = ', low_conf, 'j_conf = ', up_conf, &
   !'l_index = ', elements(el_index)%ions(current_ion)%levels(I)%l_index
   ! firstly we will check out if this is a transition which levels were really included
   found_low_conf = .FALSE.
   found_up_conf = .FALSE.
   n_levels = SIZE(elements(el_index)%ions(current_ion)%levels)
   DO J = 1, n_levels
    IF(low_conf == elements(el_index)%ions(current_ion)%levels(J)%elconf) THEN
     found_low_conf = .TRUE.
     !print*, 'found electron configuration...'
     IF(col_str >= 0) THEN
      linelist(ntransitions + 1)%upper = elements(el_index)%ions(current_ion)%levels(J)%l_index
     ELSE
      linelist(ntransitions + 1)%lower = elements(el_index)%ions(current_ion)%levels(J)%l_index
     END IF 
     !print*, 'l_index = ', elements(el_index)%ions(current_ion)%levels(J)%l_index
    END IF
    IF(up_conf == elements(el_index)%ions(current_ion)%levels(J)%elconf) THEN
     found_up_conf = .TRUE.
     !print*, 'found electron configuration...'
     IF(col_str >= 0) THEN
      linelist(ntransitions + 1)%lower = elements(el_index)%ions(current_ion)%levels(J)%l_index
     ELSE
      linelist(ntransitions + 1)%upper = elements(el_index)%ions(current_ion)%levels(J)%l_index
     END IF
     !print*, 'l_index = ', elements(el_index)%ions(current_ion)%levels(J)%l_index
    END IF
   END DO
   IF(found_low_conf .EQV. .FALSE. .OR. found_up_conf .EQV. .FALSE.) THEN
    ! this configuration will not be taken into account and we will read the next line
    print*, 'element: ', element, ' ion = ', ion_index, ' line from ', low_conf, ' to ', up_conf, &
     '  was not included...'
    CYCLE
   ELSE
    ntransitions = ntransitions + 1
   END IF
   IF((current_element /= element) .OR. (current_ion < lowerion) &
     .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC TRANSITIONS...'
   !print*, 'indexe = ', el_index, ' indexi = ', current_ion, ' ntransitions = ', ntransitions
   linelist(ntransitions)%indexe = el_index
   linelist(ntransitions)%indexi = current_ion
   linelist(ntransitions)%freq = 1.E+8*light_speed/l_freq
   linelist(ntransitions)%A_ul = abs(A)
   linelist(ntransitions)%f_ul = 1E-3 * oconstant / &
    elements(el_index)%ions(current_ion)%levels(low_level)%stat_waight *&
    (abs(linelist(ntransitions)%A_ul) / linelist(ntransitions)%freq ** 2)
   !up_conf, elements(el_index)%ions(current_ion)%levels(I)%elconf, &
   !'l_index = ', elements(el_index)%ions(current_ion)%levels(I)%l_index
   ! we have a transition between two atomic levels and we have to connect
   ! this transition with already read levels from another file, we do this
   ! connection via electron configuration in the form of string
   trans = trans + 1
   IF(trans == ntrans(I)) EXIT
  END DO
 END DO
 ! if we did now use every transition in the file we will reallocate the array
 ! linelist so it will not be so large
 IF(ntransitions < SIZE(linelist)) THEN
  ALLOCATE(pom(SIZE(linelist)))
  pom(1:new_size) = linelist(1:new_size)
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
  print*, 'no transitions found'
  STOP
 END SELECT
CLOSE(9) 

END SUBROUTINE read_transitions


