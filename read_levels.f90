! reads atomic data for the given element (from the file compose_adata.dat)
! in the format given by the parameter 
SUBROUTINE read_levels(el_index, lowerion, upperion, levels_type, filename)

  USE types

  IMPLICIT NONE    

 ! input parameters
 INTEGER                        :: element, lowerion, upperion, levels_type, el_index
 CHARACTER (LEN=20)             :: filename
 ! loop variables
 INTEGER                        :: I, J, K
 INTEGER                        :: act_lev
 ! reading from file variables
 INTEGER                        :: ios, reading_levels
 INTEGER                        :: current_element, current_ion, ions
 INTEGER                        :: n_levels, kindex, l_numb, l_index
 INTEGER                        :: indexi, act_index
 INTEGER                        :: lowering_index
 ! saves number of levels for the ions (lowerion, upperion)
 INTEGER, ALLOCATABLE           :: nlevels(:)
 INTEGER                        :: act_nlevels, nions, nline
 INTEGER                        :: cur_ion, ionindex
 CHARACTER (LEN=200)            :: line
 CHARACTER (LEN=15)              :: iconf
 CHARACTER (LEN=15)              :: junk
 DOUBLE PRECISION               :: i_pot, l_energy, ionoffset, ionstage, s_weight
 DOUBLE PRECISION, PARAMETER    :: rydberg = 13.5979996 !(eV)
 ! calculation of excitation energy (OP)
 ! current excitation energy
 DOUBLE PRECISION               :: cur_excien 
 INTEGER                        :: cur_level
 INTEGER                        :: n_ions
 INTEGER                        :: at_index
 ! basic setting of variables
 ionoffset = 0
 ions = 0
 element = elements(el_index)%atom_number
OPEN(8,status='old',FILE=filename)
 ! now we will choose the reading file using the given levels type
 SELECT CASE(levels_type)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! standard atomic data (26.03.2016)
 ! **0** the first input of the atomic data for the 3D wind code
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE(0)
  DO
   ! reading from the file filename
   ! we ignore rows starting *
   ! firstly we read the whole line
   READ(8,'(A)',IOSTAT=reading_levels) line
   IF (ios /= 0) THEN
    print*, 'the subroutine read_atomic data:'
    STOP 'ERROR: NO VALID ATOMIC DATA...'
   END IF
   IF(line(1:1) == '*') CYCLE
   ! if everything is OK, we will read from the variable line variables
   READ(line,*) current_element, current_ion, n_levels, i_pot
   print*, 'current element: ', current_element, 'current_ion: ', current_ion, &
        'number of levels: ', n_levels, 'i_pot: ', i_pot
   ! do we read the right file?
   IF((current_element /= element).OR. (current_ion < lowerion) .OR. (current_ion > upperion)) STOP 'WRONG ATOMIC DATA...'
   ! computation and allocation of important variables
   ionstage = elements(current_element)%ions(current_ion)%ion_stage
   elements(current_element)%ions(current_ion)%ion_potential = i_pot * e_v
   !elements(current_element)%ions(current_ion)%nlevels = n_levels
   ALLOCATE(elements(current_element)%ions(current_ion)%levels(n_levels))
   ! now we will read every single atomic levels
   DO J = 1, n_levels
    READ(8,*,IOSTAT=reading_levels) l_numb, l_energy, s_weight, junk
    ! did we read anything?
    IF(reading_levels /= 0) STOP 'WRONG NUMBER OF LEVELS IN THE FILE...'
    print*, 'from the level file: ', l_numb, l_energy, s_weight, junk
    ! increase the level energy l_e (multiply with e_v) by ionoffset i. e.
    ! for the ion. potential of the ground level
    elements(current_element)%ions(current_ion)%levels(J)%exci_energy = l_energy * e_v + ionoffset
    elements(current_element)%ions(current_ion)%levels(J)%stat_waight = s_weight
    print*, 'exci_energy = ', elements(current_element)%ions(current_ion)%levels(J)%exci_energy/e_v
   END DO
    ! we can increase number of read ions :-)
    ions = ions + 1
    ! are every single ions already read?
    IF (ions == (upperion - lowerion + 1)) EXIT
  END DO
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! data from the Opacity project
 ! we expect data in this form
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE(2)
 ! we now read ionization potentials for the given ions
 nions = upperion - lowerion + 1
 at_index = elements(el_index)%atom_number
! I = 0
! DO 
!  READ(8,'(A)') line
!  !print*, line
!  IF(line(1:1) == '*') CYCLE
!  READ(line,*) indexi, i_pot
!  I = I + 1
!  IF(indexi /= I) STOP 'wrong atomic data: ion indexes are not equal'
!  ! elements(el_index)%ions(I)%ion_potential  = i_pot * e_v
!  IF(I == (upperion - lowerion + 1)) EXIT
! END DO
 ! number of levels for every ion
 ALLOCATE(nlevels(nions))
 ! set up initial variables for every single index I
 DO I = 1, nions
  nlevels(I) = 0
 END DO
 ! computing number of energy levels for the given ions
 current_ion = lowerion
 DO
  READ(8,'(A)', iostat = reading_levels) line
  IF(reading_levels /= 0) EXIT
  !print*, line
  IF(line(1:1) .EQ. '*') CYCLE
  READ(line,*) junk, junk, cur_ion, junk, junk, junk, junk, junk
  ! calculation of index of the given ion of the array nlevels(:)
  act_index = cur_ion - lowerion + 2
  nlevels(act_index) = nlevels(act_index) + 1
 END DO
 !print*, 'number of levels: ', nlevels
 DO I = 1, nions
  ! ionindex = I + lowerion - 1
  indexi = at_index - I - lowerion + 3
  n_levels = nlevels(I)
  ! write(*,*) 'I = ', I, 'lowerion = ', lowerion, &
  !  ' el_index = ', el_index, ' indexi = ', indexi, ' n_levels = ', n_levels
  ALLOCATE(elements(el_index)%ions(indexi)%levels(n_levels))
 END DO
 ! allocation of the given arrays
 REWIND(8)
 ! we have to find again the data
! nline = 0
! DO
!  READ(8,'(A)') line
!  IF(line(1:1) .EQ. '*') CYCLE
!  nline = nline + 1
!  IF(nline == nions) EXIT
! END DO
 ! now we are reading atomic data for the selected ions
 DO current_ion = lowerion, upperion ! loop over ions
  ! index in the array elements%ions(indexi) ordered from
  ! neutrals to most ionized ions
  indexi = at_index - current_ion + 2
  act_index = current_ion - lowerion + 1
  act_nlevels = nlevels(act_index)
  !print*, 'n_levels = ', act_nlevels
   ! we have to calculate ionoffset
!   ionoffset = 0
!   IF(current_ion > 1) THEN
!    DO I = 1, current_ion - 1
!     ionoffset = ionoffset + elements(el_index)%ions(I)%ion_potential / e_v
!    END DO
!     !ionoffset = ionoffset + i_pot 
!   ELSE
!    ionoffset = 0
!   END IF
  ! reading the levels for the given ion
  J = 0
  DO  ! loop over atomic levels for the given ion
   READ(8,'(A)', IOSTAT = reading_levels) line
   ! print*, 'reading_levels: ', line
   IF(reading_levels /= 0) EXIT
   IF(line(1:1) .EQ. '*') CYCLE
   READ(line,*) kindex, junk, junk, junk, l_index, iconf, l_energy, s_weight
   J = J + 1
   elements(el_index)%ions(indexi)%levels(J)%exci_energy = l_energy * rydberg * e_v
   !write(*,*) 'read_levels: kindex = ', kindex, ' element = ', el_index, 'ion = ', indexi, &
   ! ' J = ', J, ' exci_energy = ', &
   ! elements(el_index)%ions(indexi)%levels(J)%exci_energy / e_v, ' l_energy = ', l_energy
   elements(el_index)%ions(indexi)%levels(J)%stat_waight = s_weight
   elements(el_index)%ions(indexi)%levels(J)%elconf = iconf
   ! we have to calculate l_index correctly: it should start at 1 for every ion
   ! this condition is not satisfied in the input files thus we have to substract
   ! the total number of levels of the lower ions from the number kindex
   IF(J == 1) lowering_index = kindex
   elements(el_index)%ions(indexi)%levels(J)%l_index = kindex - lowering_index + 1
   IF( J == act_nlevels) EXIT
  END DO ! loop over atomic levels for the given ion
  elements(el_index)%ions(indexi)%ion_potential = &
   ABS(MINVAL(elements(el_index)%ions(indexi)%levels(:)%exci_energy))
  ! write(*,*) 'read_levels: ion pot = ', MINVAL(elements(el_index)%ions(indexi)%levels(:)%exci_energy)/ e_v
  DO act_lev = 1, SIZE(elements(el_index)%ions(indexi)%levels)
   elements(el_index)%ions(indexi)%levels(act_lev)%exci_energy = &
    elements(el_index)%ions(indexi)%levels(act_lev)%exci_energy + &
    elements(el_index)%ions(indexi)%ion_potential
    ! write(*,*) 'read_levels: el = ', el_index, ' ion = ', indexi, &
    !  ' J = ', J, ' act_lev = ', act_lev, ' excie = ', &
    !  elements(el_index)%ions(indexi)%levels(act_lev)%exci_energy
  END DO
 END DO ! loop over ions
 ! recalculation of excitation energies
 ! this sbr recalculates excitation energies WRT of the ionization
 ! energy downloaded from the Opacity project 
 n_ions = SIZE(elements(el_index)%ions)
 ionoffset = 0.D0
 DO I = 1, n_ions
  n_levels = SIZE(elements(el_index)%ions(I)%levels)
  IF(I > 1) ionoffset = ionoffset + elements(el_index)%ions(I - 1)%ion_potential
  ! write(*,*) 'read_levels: ion pot = ', elements(el_index)%ions(I)%ion_potential/ e_v
  DO cur_level = 1, n_levels
   cur_excien = elements(el_index)%ions(I)%levels(cur_level)%exci_energy
   ! write(*,*) 'read_levels: cur_excien = ', cur_excien,&
   !  ' ee = ', elements(el_index)%ions(I)%levels(cur_level)%exci_energy
   elements(el_index)%ions(I)%levels(cur_level)%exci_energy = ionoffset + cur_excien
   ! write(*,*) 'read_levels: el_index = ', el_index, 'cur_ion = ', I, &
   !  ' cur_level = ', cur_level, 'ionoffset = ', ionoffset / e_v, 'eenergy = ', &
   !  elements(el_index)%ions(I)%levels(cur_level)%exci_energy / e_v
  END DO
 END DO
  
   ! if everything is OK, we will read from the variable line variables
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! **DEFAULT** default case
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE DEFAULT
  print*, 'choice has not been found'
  print*, 'reading only hydrogen data inluded in the code...'
  STOP 'with no valid atomic data...'
 END SELECT
CLOSE(8)
END SUBROUTINE read_levels
