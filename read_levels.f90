! reads atomic data for the given element (from the file compose_adata.dat)
! in the format given by the parameter 
SUBROUTINE read_levels(el_index, lowerion, upperion, levels_type, filename)

  USE types

  IMPLICIT NONE    

 ! input parameters
 INTEGER                        :: element, lowerion, upperion, levels_type, el_index
 INTEGER                        :: levelindex
 CHARACTER (LEN=20)             :: filename
 ! loop variables
 INTEGER                        :: I, J
 INTEGER                        :: act_lev
 ! reading from file variables
 INTEGER                        :: reading_levels
 INTEGER                        :: current_ion, ions
 INTEGER                        :: n_levels
 INTEGER                        :: indexi, act_index
 ! saves number of levels for the ions (lowerion, upperion)
 INTEGER, ALLOCATABLE           :: nlevels(:)
 INTEGER                        :: act_nlevels, nions
 INTEGER                        :: cur_ion, vsplit
 CHARACTER (LEN=200)            :: line
 CHARACTER (LEN=16)              :: iconf
 CHARACTER (LEN=16)              :: junk
 INTEGER                        :: jint
 REAL                           :: jreal
 DOUBLE PRECISION               :: l_energy, ionoffset, s_weight
 DOUBLE PRECISION, PARAMETER    :: rydberg = 13.5979996 !(eV)
 ! calculation of excitation energy (OP)
 ! current excitation energy
 DOUBLE PRECISION               :: cur_excien 
 INTEGER                        :: cur_level
 INTEGER                        :: n_ions
 INTEGER                        :: at_index
 ! basic setting of variables
! format of Opacity project file
991 format(I7, I3, I3, I5, I4, Tr1, A16, e12.5, f5.1)
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
  STOP 'reading_levels: this level sources is not possible'
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
!  !write(99,*) line
!  IF(line(1:1) == '*') CYCLE
!  READ(line,*) indexi, i_pot
!  I = I + 1
!  IF(indexi /= I) STOP 'wrong atomic data: ion indexes are not equal'
!  ! elements(el_index)%ions(I)%ion_potential  = i_pot * e_v
!  IF(I == (upperion - lowerion + 1)) EXIT
! END DO
 ! number of levels for every ion
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! it must be changed because data for an atom can be stored in severe files
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ALLOCATE(nlevels(nions))
 ! set up initial variables for every single index I
 DO I = 1, nions
  nlevels(I) = 0
 END DO
 ! computing number of energy levels for the given ions
 current_ion = lowerion
 DO
  !
  ! calculation of number of levels for the given ion
  !
  ! READ(8,'(A)', iostat = reading_levels) line
  ! write(*,*) 'read_levels: line = ', line
  ! BACKSPACE(8)
  READ(8 ,991, iostat = reading_levels) jint, jint, cur_ion, jint, jint, junk, l_energy, jreal
  ! write(*,*) 'reading_levels = cur_ion = ', cur_ion, ' l_energy = ', l_energy
  IF(reading_levels /= 0) EXIT
  ! we ignore energy levels in the continuum
  IF(l_energy > 0.D0) CYCLE
  ! calculation of index of the given ion of the array nlevels(:)
  act_index = cur_ion - lowerion + 2
  nlevels(act_index) = nlevels(act_index) + 1
 END DO
 ! write(*,*) 'number of levels: ', nlevels(act_index)
 DO I = 1, nions
  ! ionindex = I + lowerion - 1
  indexi = at_index - I - lowerion + 3
  n_levels = nlevels(I)
  ! write(99,*) 'I = ', I, 'lowerion = ', lowerion, &
  !  ' el_index = ', el_index, ' indexi = ', indexi
  ! write(*,*) 'read_levels: ion = ', I, ' n_levels = ', n_levels
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
  !write(99,*) 'n_levels = ', act_nlevels
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
   ! READ(8,'(A)', IOSTAT = reading_levels) line
   ! write(99,*) 'reading_levels: ', line
   ! finally read the atomic data and save it into the variables
   READ(8 ,991, iostat = reading_levels) levelindex, jint, jint, vsplit, jint, iconf, l_energy, s_weight
   IF(reading_levels /= 0) EXIT
   ! IF(line(1:1) .EQ. '*') CYCLE
   ! READ(line,*) levelindex, junk, junk, junk, l_index, iconf, l_energy, s_weight
   IF(l_energy > 0.D0) CYCLE
   J = J + 1
   elements(el_index)%ions(indexi)%levels(J)%exci_energy = l_energy * rydberg * e_v
   elements(el_index)%ions(indexi)%levels(J)%levelindex = levelindex
   ! write(99,*) 'read_levels: kindex = ', kindex, ' element = ', el_index, 'ion = ', indexi, &
   !  ' J = ', J, ' exci_energy = ', &
   ! elements(el_index)%ions(indexi)%levels(J)%exci_energy / e_v, ' l_energy = ', l_energy
   ! write(*,*) 'read_levels: element = ', el_index, 'ion = ', indexi, &
   !  ' J = ', J, ' exci_energy = ', &
   !  elements(el_index)%ions(indexi)%levels(J)%exci_energy / e_v, ' l_energy = ', l_energy
   elements(el_index)%ions(indexi)%levels(J)%stat_waight = s_weight
   elements(el_index)%ions(indexi)%levels(J)%elconf = iconf
   elements(el_index)%ions(indexi)%levels(J)%vsplit = vsplit
   ! we have to calculate l_index correctly: it should start at 1 for every ion
   ! this condition is not satisfied in the input files thus we have to substract
   ! the total number of levels of the lower ions from the number kindex
   ! elements(el_index)%ions(indexi)%levels(J)%l_index = J
   IF( J == act_nlevels) EXIT
  END DO ! loop over atomic levels for the given ion
  elements(el_index)%ions(indexi)%ion_potential = &
   ABS(MINVAL(elements(el_index)%ions(indexi)%levels(:)%exci_energy))
  ! write(*,*) 'read_levels: ion pot = ', MINVAL(elements(el_index)%ions(indexi)%levels(:)%exci_energy)/ e_v
  DO act_lev = 1, SIZE(elements(el_index)%ions(indexi)%levels)
   elements(el_index)%ions(indexi)%levels(act_lev)%exci_energy = &
   elements(el_index)%ions(indexi)%levels(act_lev)%exci_energy + &
    elements(el_index)%ions(indexi)%ion_potential
    ! write(99,*) 'read_levels: el = ', el_index, ' ion = ', indexi, &
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
   ! IF(el_index == 3) write(*,*) 'read_levels: el_index = ', el_index, 'cur_ion = ', I, &
   !  ' cur_level = ', cur_level, 'ionoffset = ', ionoffset / e_v, 'eenergy = ', &
   !  elements(el_index)%ions(I)%levels(cur_level)%exci_energy / e_v
   ! write(*,*) 'read_levels: exci_energy = ', elements(el_index)%ions(I)%levels(cur_level)%exci_energy
  END DO
 END DO
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! data from the Tardis
 ! we expect data in this form
 ! atomic_number ion_number level_number  energy   g  metastable
 ! we won't use the metastable variable
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE(8)
  nions = upperion - lowerion + 1
  at_index = elements(el_index)%atom_number
  ALLOCATE(nlevels(nions))
  ! set up initial variables for every single index I
  DO I = 1, nions
   nlevels(I) = 0
  END DO
  ! computing number of energy levels for the given ions
  current_ion = lowerion
  ! reading the first line
  READ(8 ,*, iostat = reading_levels) jint, jint, iconf, l_energy, s_weight, junk
  DO
   !
   ! calculation of number of levels for the given ion
   !
   ! READ(8,'(A)', iostat = reading_levels) line
   ! write(*,*) 'read_levels: line = ', line
   ! BACKSPACE(8)
   READ(8 ,*, iostat = reading_levels) jint, cur_ion, jint, jreal, jint, junk
   IF(reading_levels /= 0) EXIT
   act_index = cur_ion + 1
   nlevels(act_index) = nlevels(act_index) + 1
  END DO
  DO I = 1, nions
   ! ionindex = I + lowerion - 1
   indexi = I + lowerion - 1
   n_levels = nlevels(I)
   ! write(99,*) 'I = ', I, 'lowerion = ', lowerion, &
   !  ' el_index = ', el_index, ' indexi = ', indexi
   ! write(*,*) 'read_levels: ion = ', I, ' n_levels = ', n_levels
   ALLOCATE(elements(el_index)%ions(indexi)%levels(n_levels))
  END DO
  ! allocation of the given arrays
  REWIND(8)
  ! we have to find again the data
!  nline = 0
!  DO
!   READ(8,'(A)') line
!   IF(line(1:1) .EQ. '*') CYCLE
!   nline = nline + 1
!   IF(nline == nions) EXIT
!  END DO
  ! now we are reading atomic data for the selected ions
  ! reading the first line
  READ(8 ,*, iostat = reading_levels) jint, jint, iconf, l_energy, s_weight, junk
  DO current_ion = lowerion, upperion ! loop over ions
   ! index in the array elements%ions(indexi) ordered from
   ! neutrals to most ionized ions
   indexi = current_ion
   act_index = current_ion - lowerion + 1
   act_nlevels = nlevels(act_index)
   !write(99,*) 'n_levels = ', act_nlevels
    ! we have to calculate ionoffset
!    ionoffset = 0
!    IF(current_ion > 1) THEN
!     DO I = 1, current_ion - 1
!      ionoffset = ionoffset + elements(el_index)%ions(I)%ion_potential / e_v
!     END DO
!      !ionoffset = ionoffset + i_pot 
!    ELSE
!     ionoffset = 0
!    END IF
   ! reading the levels for the given ion
   J = 0
   DO  ! loop over atomic levels for the given ion
    ! READ(8,'(A)', IOSTAT = reading_levels) line
    ! write(99,*) 'reading_levels: ', line
    ! finally read the atomic data and save it into the variables
    READ(8 ,*, iostat = reading_levels) jint, jint, iconf, l_energy, s_weight, junk
    IF(reading_levels /= 0) EXIT
    ! IF(line(1:1) .EQ. '*') CYCLE
    J = J + 1
    ! write(*,*) 'read_levels: el_index = ', el_index, ' J = ', J
    ! write(*,*) 'read_levels: iconf = ', iconf, ' l_energy = ', l_energy
    ! write(*,*) 'read_levels: indexi = ', indexi, ' act_index = ', act_index, 'act_nlevels = ', act_nlevels
    elements(el_index)%ions(indexi)%levels(J)%exci_energy = l_energy
    elements(el_index)%ions(indexi)%levels(J)%elconf = iconf
    ! write(99,*) 'read_levels: kindex = ', kindex, ' element = ', el_index, 'ion = ', indexi, &
    !  ' J = ', J, ' exci_energy = ', &
    ! elements(el_index)%ions(indexi)%levels(J)%exci_energy / e_v, ' l_energy = ', l_energy
    ! write(*,*) 'read_levels: element = ', el_index, 'ion = ', indexi, &
    !  ' J = ', J, ' exci_energy = ', &
    !  elements(el_index)%ions(indexi)%levels(J)%exci_energy / e_v, ' l_energy = ', l_energy
    elements(el_index)%ions(indexi)%levels(J)%stat_waight = s_weight
    ! we have to calculate l_index correctly: it should start at 1 for every ion
    ! this condition is not satisfied in the input files thus we have to substract
    ! the total number of levels of the lower ions from the number kindex
    ! elements(el_index)%ions(indexi)%levels(J)%l_index = J
    IF( J == act_nlevels) EXIT
   END DO ! loop over atomic levels for the given ion
   elements(el_index)%ions(indexi)%ion_potential = &
    ABS(MAXVAL(elements(el_index)%ions(indexi)%levels(:)%exci_energy))
   ! write(*,*) 'read_levels: ion pot = ', MINVAL(elements(el_index)%ions(indexi)%levels(:)%exci_energy)/ e_v
  END DO ! loop over ions
  ! recalculation of excitation energies
  ! this loop recalculates energies so every energy level will have energy equal to
  ! E = exci_energy + ion_energy
  n_ions = SIZE(elements(el_index)%ions)
  ionoffset = 0.D0
  DO I = 1, n_ions
   n_levels = SIZE(elements(el_index)%ions(I)%levels)
   IF(I > 1) ionoffset = ionoffset + elements(el_index)%ions(I - 1)%ion_potential
   ! write(*,*) 'read_levels: ion pot = ', elements(el_index)%ions(I)%ion_potential/ e_v
   DO cur_level = 1, n_levels
    cur_excien = elements(el_index)%ions(I)%levels(cur_level)%exci_energy
    elements(el_index)%ions(I)%levels(cur_level)%exci_energy = ionoffset + cur_excien
    ! write(*,*) 'read_levels: ', el_index, I, cur_level, 'ionoffset = ', ionoffset,&
    !  ' ee = ', elements(el_index)%ions(I)%levels(cur_level)%exci_energy
   END DO
  END DO
   ! if everything is OK, we will read from the variable line variables
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! **DEFAULT** default case
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 CASE DEFAULT
  write(99,*) 'choice has not been found'
  write(99,*) 'reading only hydrogen data inluded in the code...'
  STOP 'with no valid atomic data...'
 END SELECT
CLOSE(8)
END SUBROUTINE read_levels
