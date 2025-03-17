! Read the input file and saves composition, atomic levels and atomic transitions
! name of an input file can be found in the input file under a variable 'inputmodelFile'
! and it is saved into a variable 'inputmodelFile'
! 
! structure of the file
! 
! **composition**
! proton number, abundance, lower ion, upper ion, mass
! 
! **levels**
! proton number, lower ion, upper ion, type of levels, level file
! 
! **transitions**
! proton number, lower ion, upper ion, type of transitions, transitions file, photcross type, n of cross section, photcross file
!
! INPUT: NONE
! OUTPUT: NONE
! 
SUBROUTINE read_composition()


USE types
USE constants

IMPLICIT NONE    

INTEGER                            :: ind_I, ind_J
INTEGER                            :: element_index, Z, lowerion, upperion
INTEGER                            :: photn
INTEGER                            :: levels_type, transition_type
INTEGER                            :: atom_number
DOUBLE PRECISION                   :: abundance
INTEGER                            :: current_ion, nions,ios
CHARACTER (filename_length)                     :: filename, photfile
CHARACTER (LEN=filename_length)                    :: line
DOUBLE PRECISION                   :: mass, tot_abundance
! photon cross section data type
INTEGER                            :: phcs_type
LOGICAL                            :: fileExists

! initialization the toal number of photoionization cross section
n_photcrossect = 0
inquire( file=inputcomposition, exist=fileExists )
if(.NOT. fileExists) THEN
 write(*,*) "file: ", inputcomposition, "does not exists"
 STOP
end if
OPEN (UNIT=7, FILE=inputcomposition)
 ! computes number of lines in the input file
 ! number of rows is equal to 0
 ! for this time it will calculate number of rows
 n_elements = 0
 tot_abundance = 0.D0
 DO
  READ(7,*,IOSTAT=ios) line
   IF (ios /= 0) EXIT
   IF ( line == '**levels**' ) EXIT
   IF ( line(1:1) == '*') CYCLE
  n_elements = n_elements + 1
 END DO
 ! Allocate the memory to the elements(n_elements)
 write(99,*) 'number of elements: ', n_elements
 ALLOCATE (elements(n_elements)) 
 write(99,*) 'dimension of elements = ', SIZE(elements)
 REWIND(7)
 ! Loop over all rows involved i.e. read all other lines in the compose_adata.dat
 ! and assine these values to the elements(I)%... and elements(I)%ions(J)%...
  ind_I=1
  write(99,*) 'reading chemical composition'
 DO 
  READ(7,'(A)',iostat=ios) line
  IF (ios /= 0) EXIT
   IF ( TRIM(line) == '**levels**' ) THEN
      !write(99,*) 'read_composition: we have found the string **levels**...'
      EXIT
   END IF
   IF ( line(1:1) == '*') CYCLE
   READ(line,*) Z, abundance, lowerion, upperion, mass
   write(99,*) Z, abundance, lowerion, upperion, mass
   elements(ind_I)%atom_number = Z
   elements(ind_I)%atom_mass = mass * const_mp_g
   ! Number of ions 
   nions =  upperion - lowerion + 1
   !write(99,*) 'Z = ', Z, ' nions = ', nions
   elements(ind_I)%nions = nions
   elements(ind_I)%abundance = abundance
   tot_abundance = tot_abundance + abundance
   ! Assine lowerion to the current ion which we will use to caunt number of ions
   ! This is important because we can play only with 3 and 4 ion.stage of some element
   current_ion = lowerion
   !write(99,*) 'current ion =', current_ion
   !write(99,*) 'upperion = ', upperion, 'lowerion = ', lowerion, 'upperion - lowerion + 1', nions
   ! Allocate the memory to the elements(I)%ions(nions)
   ALLOCATE (elements(ind_I)%ions(nions))
   ! Loop over all ions of given chem. element
   DO ind_J = lowerion, upperion
    elements(ind_I)%ions(ind_J)%ion_stage = current_ion
    current_ion = lowerion + 1
   END DO
   ind_I = ind_I + 1
 END DO
 ! norma of abundances
 DO ind_I = 1, n_elements
  elements(ind_I)%abundance = elements(ind_I)%abundance / tot_abundance
 END DO
 ! now reading atomic levels
 write(99,*) 'reading atomic levels'
 DO 
    READ(7,'(A)',iostat=ios) line
    ! write(99,*) line
  IF (ios /= 0) EXIT
  IF ( TRIM(line) == '**transitions**' ) THEN
      EXIT
  END IF
  IF ( INDEX(line, '*') /= 0) CYCLE
  ! write(*,*) 'read_composition: read line = ', line
  READ(line,*) atom_number, lowerion, upperion, levels_type, filename
  ! write(*,*) 'read_composition: filename = ', filename
  CALL find_element_index(atom_number, element_index)
  ! the most important is to read the file
  ! write(*,*) 'read_composition: element_index = ', element_index
  CALL read_levels(element_index, lowerion,upperion,levels_type,filename)
 END DO
! now reading atomic transitions 
! ! now we can read informations if the files
 write(99,*) 'reading atomic transitions'
 IF(.NOT. oneline) THEN
  DO 
   READ(7,'(A)',iostat=ios) line
   IF (ios /= 0) EXIT
   IF ( INDEX(line, '*') /= 0) CYCLE
   READ(line,*) atom_number, lowerion, upperion, transition_type, filename, &
         phcs_type, photn, photfile
   !write(99,*) 'calling subroutine read_transitions...'
   CALL find_element_index(atom_number, element_index)
   CALL read_transitions(element_index, lowerion, upperion, transition_type, filename)
   IF(photn == 0) THEN
    write(99,*) 'no valid data for potoionization cross sections'
   ELSE
    write(99,*) 'calling subroutine read_photcs... for element index = ', element_index
    CALL read_photcs(phcs_type, element_index, photn, photfile)
   END IF
  END DO
 ELSE IF(oneline) THEN
  CALL read_singleline()
 END IF
 
! OPEN(20,status='new',FILE='linelist.dat')
!  DO ind_I = 1, ntransitions
!   write(20,*) elements(linelist(ind_I)%indexe)%atom_number, linelist(ind_I)%indexi, linelist(ind_I)%freq, linelist(ind_I)%f_ul
!  END DO
! CLOSE(20)
 ! sorting the linelist ray 
 write(99,*) 'read_composition: dim of linelist = ', SIZE(linelist)
 ! STOP 'read_composition, testing'
 write(99,*) 'sorting the linelist'
 CALL sorting_new(ntransitions, linelist)
  !write(99,*) 'ntransitions = ', ntransitions
 !!! Only for testing
! PRINT*, 'testing'
! DO ind_I = 1, n_elements
!  element_index = I
!  nions = SIZE(elements(ind_I)%ions)
!  Z = elements(ind_I)%atom_number
!  lowerion = elements(ind_I)%ions(1)%ion_stage
!  upperion = elements(ind_I)%ions(nions)%ion_stage
!  PRINT*, element_index, Z, lowerion, upperion
! END DO
! DO ind_I = 1, n_elements
!  n_ions = SIZE(elements(ind_I)%ions)
!  DO J = 1, n_ions
!   n_levels = SIZE(elements(ind_I)%ions(J)%levels)
!   DO K = 1, n_levels
!    n_points = SIZE(elements(ind_I)%ions(J)%levels(K)%photcros)
!    write(99,*) 'element = ', I, ' ion = ', J, ' level = ', K, ' n_points = ', n_points
!   END DO
!  END DO
! END DO
! OPEN(20,status='new',FILE='linelist_ordered.dat')
  DO ind_I = 1, ntransitions
   linelist(ind_I)%n_deexc = 0
   linelist(ind_I)%n_exc = 0
   ! write(*,*) 'read_composition: line = ', I, ' lower = ', linelist(ind_I)%lower, ' upper = ', linelist(ind_I)%upper
!   write(20,*) elements(linelist(ind_I)%indexe)%atom_number, linelist(ind_I)%indexi, linelist(ind_I)%freq, linelist(ind_I)%f_ul
  END DO
! CLOSE(20)
! STOP 'read_composition: testing'

END SUBROUTINE read_composition
