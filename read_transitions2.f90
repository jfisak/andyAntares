SUBROUTINE read_transitions2(element, lowerion, upperion, transition_type, filename)


 USE types

 IMPLICIT NONE    

 ! input data
 INTEGER                        :: element, lowerion, upperion, transition_type
 INTEGER                        :: n_ions, junk
 CHARACTER (LEN=20)             :: filename
 ! used constants
 DOUBLE PRECISION               :: oconstant
 ! loop variables
 INTEGER                        :: I
 ! reading from the file
 CHARACTER (LEN=200)            :: line
 INTEGER                        :: ios, linereading, n_levels
 INTEGER                        :: n_transitions, n_line
 INTEGER                        :: current_element, current_ion
 INTEGER                        :: kindex, low_level, up_level
 CHARACTER (LEN=4)              :: low_conf, up_conf
 DOUBLE PRECISION               :: A, col_str, l_freq


! calculate the constant for the oscilator strength calculation
oconstant = (me_g * light_speed ** 3)/(8.D0 * pi ** 2 * e_charge**2)
n_ions = 0

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
  n_ions = elements(current_element)%nions
  print*, 'number of ions n_ions = ', n_ions
   DO I = 1, n_transitions
    READ(9,'(A)',IOSTAT=linereading) line
     IF (linereading /= 0) THEN
     print*, 'the subroutine read_transitions2:'
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
    linelist(n_line)%freq = (elements(current_element)%ions(current_ion)%levels(up_level)%exci_energy -&
                               elements(current_element)%ions(current_ion)%levels(low_level)%exci_energy) / h
    linelist(n_line)%A_ul = A
    linelist(n_line)%f_ul = oconstant * (elements(current_element)%ions(current_ion)%levels(up_level)%stat_waight / &
                       elements(current_element)%ions(current_ion)%levels(low_level)%stat_waight) *                 &
                       (linelist(n_line)%A_ul / linelist(n_line)%freq ** 2)
    print*, n_line, ': ', 'indexe: ', linelist(n_line)%indexe, linelist(n_line)%indexi, linelist(n_line)%lower, &
                       linelist(n_line)%upper, linelist(n_line)%freq, linelist(n_line)%A_ul, linelist(n_line)%f_ul
   END DO
   n_ions = n_ions + 1
   IF(n_ions == (upperion - lowerion + 1)) EXIT
 END DO
 ! if we do not have the ray linelist allocated, we have
 ! to calculate number of possible transitions
 ELSE
 DO   
  READ(9,'(A)',IOSTAT=ios) line
  print*, line
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
! we expect data in this form
! **iont**
! Z    iont    n_trans
! **lines**
! i    iLV     jLV     iCONF   jCONF   gF      gA      WL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE(2)
  ! only if the variable linelis is allocated
  IF(ALLOCATED(linelist)) THEN
   ! do 0
   DO 
    ! do 1
    ! looking for the flag **iont**
    DO
     READ(9,'(A)',IOSTAT=ios) line
     print*, line
     IF (ios /= 0) EXIT
     IF(TRIM(line) == '**iont**') EXIT
    ! end do 1
    END DO
    ! do 5
    DO
      READ(9,'(A)',IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF ( INDEX(line, '*') /= 0) CYCLE
      READ(line,*) current_element, current_ion, n_transitions
      IF(n_transitions == 0) CYCLE
      ! do we have the right file?
      IF((current_element /= element) .OR. (current_ion < lowerion) &
        .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC TRANSITIONS...'
     ! do 2
     DO
      READ(9,'(A)',IOSTAT=ios) line
      IF (ios /= 0) STOP 'NO VALID DATA IN THE FILE...'
      IF(TRIM(line) == '**lines**') EXIT
     ! end do 2
     END DO
     ! do 3
     DO
      READ(9,'(A)',IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF ( INDEX(line, '*') /= 0) CYCLE
      READ(line,*) kindex, low_level, up_level, low_conf, up_conf, col_str, A, l_freq
      n_line = ntransitions + 1
      ntransitions = n_line
      n_levels = elements(current_element)%ions(current_ion)%nlevels
      ! do 4
      ! will find lower and upper index for every transition
      DO I = 1,n_levels
       IF(low_conf == elements(current_element)%ions(current_ion)%levels(I)%elconf) &
          linelist(n_line)%lower = elements(current_element)%ions(current_ion)%levels(I)%l_index
       IF(up_conf == elements(current_element)%ions(current_ion)%levels(I)%elconf) &
          linelist(n_line)%upper = elements(current_element)%ions(current_ion)%levels(I)%l_index
      ! end do 4
      END DO
      ! 
      linelist(n_line)%freq = 1.E-10*light_speed/l_freq
      linelist(n_line)%f_ul = col_str
      linelist(n_line)%A_ul = A
     ! end do 3
     END DO
     IF(ios /= 0) EXIT
     n_ions = n_ions + 1
     IF(n_ions == (upperion - lowerion + 1)) EXIT
    ! end do 5
    END DO   
    ! end do 0
    END DO
  ! if we do not have the ray linelist allocated, we have
  ! to calculate number of possible transitions
  ELSE
    DO
     READ(9,'(A)',IOSTAT=ios) line
     IF (ios /= 0) STOP 'NO VALID DATA IN THE FILE...'
     IF(TRIM(line) == '**iont**') EXIT
    END DO
    ! we will compute number of lines in the opened file
    DO
     READ(9,'(A)',IOSTAT=ios) line
     print*, line
     IF (ios /= 0) EXIT
     IF ( INDEX(line, '*') /= 0) CYCLE
     READ(line,*) junk, junk, n_transitions
     ntransitions = ntransitions + n_transitions
     IF(n_transitions == 0) CYCLE
     ! loop which will find another label **iont**
     DO 
      READ(9,'(A)',IOSTAT=ios) line
      IF (ios /= 0) EXIT
      IF(TRIM(line) == '**iont**') EXIT
     END DO
    END DO
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
END SUBROUTINE read_transitions2


