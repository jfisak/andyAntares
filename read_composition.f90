SUBROUTINE read_composition()

! Read the file compose_adata.dat and remember
! elements(n_elements)
! elements(I)%atom_number, I=1,n_elements
! elements(I)%nions
! elements(I)%ions(nions)
! elements(I)%ions(J)%ion_stage, J=1,nions

  USE types

  IMPLICIT NONE    

  INTEGER                            :: I, J, flag
  INTEGER                            :: element_index, Z, lowerion, upperion
  INTEGER                            :: photn
  INTEGER                            :: levels_type, transition_type
  INTEGER                            :: atom_number
  DOUBLE PRECISION                   :: abundance
  INTEGER                            :: current_element,current_ion, nions,ios, NR
  CHARACTER (20)                     :: filename, photfile
  CHARACTER (LEN=200)                    :: line
  CHARACTER (1)                      :: junk
  DOUBLE PRECISION                   :: mass, tot_abundance
  ! photon cross section data type
  INTEGER                            :: phcs_type
  INTEGER                            :: K, n_levels, n_ions, n_points

OPEN (UNIT=7, FILE='compose_adata.dat')
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
 print*, 'number of elements: ', n_elements
 ALLOCATE (elements(n_elements)) 

 REWIND(7)
 ! Loop over all rows involved i.e. read all other lines in the compose_adata.dat
 ! and assine these values to the elements(I)%... and elements(I)%ions(J)%...
  I=1
 DO 
    READ(7,'(A)',iostat=ios) line
   IF (ios /= 0) EXIT
    IF ( TRIM(line) == '**levels**' ) THEN
       !print*, 'read_composition: we have found the string **levels**...'
       EXIT
    END IF
    IF ( line(1:1) == '*') CYCLE
    READ(line,*) Z, abundance, lowerion, upperion, mass
    PRINT*, Z, abundance, lowerion, upperion, mass
    elements(I)%atom_number = Z
    elements(I)%atom_mass = mass * mp_g
    ! Number of ions 
    nions =  upperion - lowerion + 1
    !print*, 'Z = ', Z, ' nions = ', nions
    elements(I)%nions = nions
    elements(I)%abundance = abundance
    tot_abundance = tot_abundance + abundance
    ! Assine lowerion to the current ion which we will use to caunt number of ions
    ! This is important because we can play only with 3 and 4 ion.stage of some element
    current_ion = lowerion
    !print*, 'current ion =', current_ion
    !print*, 'upperion = ', upperion, 'lowerion = ', lowerion, 'upperion - lowerion + 1', nions
    ! Allocate the memory to the elements(I)%ions(nions)
    ALLOCATE (elements(I)%ions(nions))
    if(ALLOCATED(elements(I)%ions)) print*, 'allocated: elements(', I, ')%ions...', nions
    ! Loop over all ions of given chem. element
    DO J = lowerion, upperion
       elements(I)%ions(J)%ion_stage = current_ion
       current_ion = lowerion + 1
    END DO
    I=I+1
 END DO
 ! norma of abundances
 DO I = 1, n_elements
  elements(I)%abundance = elements(I)%abundance / tot_abundance
 END DO
 ! now reading atomic levels
 DO 
    READ(7,'(A)',iostat=ios) line
    ! print*, line
  IF (ios /= 0) EXIT
  IF ( TRIM(line) == '**transitions**' ) THEN
      !print*, 'we have found **transitions**...'
      EXIT
  END IF
  IF ( INDEX(line, '*') /= 0) CYCLE
  READ(line,*) atom_number, lowerion, upperion, levels_type, filename
  CALL find_element_index(atom_number, element_index)
  ! the most important is to read the file
  CALL read_levels(element_index, lowerion,upperion,levels_type,filename)
 END DO
! now reading atomic transitions 
! ntransitions = 0
 ! now we have to compute number of possible transitions
 ! as always we skip over rows starting with '*'
! DO 
!    READ(7,'(A)',iostat=ios) line
!  IF (ios /= 0) EXIT
!    print*, line
!  IF ( INDEX(line, '*') /= 0) CYCLE
!  READ(line,*) element_index, current_element, lowerion, upperion, transition_type, filename
!  ! the most important is to read the file
!  CALL read_transitions(element_index, current_element,lowerion,upperion,transition_type,filename)
! END DO
! ! now we can allocate ray linelist
! ALLOCATE(linelist(ntransitions))
!  print*, 'ntransitions = ', ntransitions
!  ! temporary characteristics of number of saved lines
!  !linelist(ntransitions)%indexe = 0
!  ntransitions = 0
!  print*, 'ntransitions = ', ntransitions
!  IF(ALLOCATED(linelist)) print*, 'linelist was allocated succefully...'
!  ! this is used for information about fullfiled number of elements in array linelist
!  REWIND(7)
!  ! we have to find a flag **transitions**
! DO
!    READ(7,'(A)',iostat=ios) line
!    IF ( TRIM(line) == '**transitions**' ) EXIT
! END DO
! ! now we can read informations if the files
 DO 
  READ(7,'(A)',iostat=ios) line
  IF (ios /= 0) EXIT
  IF ( INDEX(line, '*') /= 0) CYCLE
  READ(line,*) atom_number, lowerion, upperion, transition_type, filename, &
        phcs_type, photn, photfile
  print*, 'calling subroutine read_transitions...'
  CALL find_element_index(atom_number, element_index)
  CALL read_transitions(element_index, lowerion, upperion, transition_type, filename)
  IF(photn == 0) THEN
   print*, 'no valid data for potoionization cross sections'
  ELSE
   print*, 'calling subroutine read_photcs...'
   CALL read_photcs(phcs_type, element_index, photn, photfile)
  END IF
 END DO
! OPEN(20,status='new',FILE='linelist.dat')
!  DO I = 1, ntransitions
!   write(20,*) elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, linelist(I)%freq, linelist(I)%f_ul
!  END DO
! CLOSE(20)
 ! sorting the linelist ray 
 CALL sorting_new(ntransitions, linelist)
  !print*, 'ntransitions = ', ntransitions
 !!! Only for testing
! PRINT*, 'testing'
! DO I = 1, n_elements
!  element_index = I
!  nions = SIZE(elements(I)%ions)
!  Z = elements(I)%atom_number
!  lowerion = elements(I)%ions(1)%ion_stage
!  upperion = elements(I)%ions(nions)%ion_stage
!  PRINT*, element_index, Z, lowerion, upperion
! END DO
! DO I = 1, n_elements
!  n_ions = SIZE(elements(I)%ions)
!  DO J = 1, n_ions
!   n_levels = SIZE(elements(I)%ions(J)%levels)
!   DO K = 1, n_levels
!    n_points = SIZE(elements(I)%ions(J)%levels(K)%photcros)
!    print*, 'element = ', I, ' ion = ', J, ' level = ', K, ' n_points = ', n_points
!   END DO
!  END DO
! END DO
! OPEN(20,status='new',FILE='linelist_ordered.dat')
!  DO I = 1, ntransitions
!   write(20,*) elements(linelist(I)%indexe)%atom_number, linelist(I)%indexi, linelist(I)%freq, linelist(I)%f_ul
!  END DO
! CLOSE(20)

END SUBROUTINE read_composition
