SUBROUTINE read_composition()

! Read the file compose_adata.dat and remember
! elements(n_elements)
! elements(I)%atom_number, I=1,n_elements
! elements(I)%nions
! elements(I)%ions(nions)
! elements(I)%ions(J)%ion_stage, J=1,nions

  USE types

  IMPLICIT NONE    

  INTEGER            :: I, J, flag
  INTEGER            :: element_index, Z, lowerion, upperion
  INTEGER            :: current_ion, nions,ios, NR
  INTEGER, PARAMETER :: maxelements = 180
  CHARACTER (20)     :: elementfile
  CHARACTER (1)      :: junk
  DOUBLE PRECISION   :: mass
  INTEGER, ALLOCATABLE :: elindexes(:)

  OPEN (UNIT=7, FILE='compose_adata.dat')
  ! computes number of lines in the input file
  ! number of rows is equal to 0
  ! for this time it will calculate number of rows
  NR = 0
  DO I=1,maxelements
   READ(7,*,IOSTAT=ios) junk, junk, junk, junk, junk, junk
    IF (ios /= 0) EXIT
    IF (I == maxelements) THEN
     print*, 'Error: Maximum number of records exceeded...'
     print*, 'Exiting program now...'
     STOP
    END IF
   NR = NR + 1
  END DO
  REWIND(7)
  print*, 'number of rows is equal to ', NR
  ALLOCATE(elindexes(NR))
  ! number of elements calculation
  DO I = 1, NR
     READ(7,*) element_index, Z, lowerion, upperion, mass, elementfile
  END DO
 IF (NR .EQ. 1) THEN
  n_elements = 1
 ELSE
  n_elements=0
  DO I = 2,NR
   flag = 0
   DO J = 1,I
    IF(elindexes(J) .EQ. elindexes(I)) flag = 1
   END DO
    IF(flag .EQ. 0) n_elements = n_elements + 1
  END DO
 END IF
  ! Allocate the memory to the elements(n_elements)
  print*, 'number of elements: ', n_elements
  ALLOCATE (elements(n_elements)) 

  REWIND(7)
  ! Loop over all rows involved i.e. read all other lines in the compose_adata.dat
  ! and assine these values to the elements(I)%... and elements(I)%ions(J)%...
  DO I = 1, NR
     READ(7,*) element_index, Z, lowerion, upperion, mass, elementfile
     PRINT*, element_index, Z, lowerion, upperion, mass
     elements(I)%atom_number = Z
     elements(I)%atom_mass = mass * mp_g
     elements(I)%levelfile = elementfile
     write(*,*) elementfile
     ! Number of ions 
     nions =  upperion - lowerion + 1
     elements(I)%nions = nions
     ! Assine lowerion to the current ion which we will use to caunt number of ions
     ! This is important because we can play only with 3 and 4 ion.stage of some element
     current_ion = lowerion
     print*, 'current ion =', current_ion  
     ! Allocate the memory to the elements(I)%ions(nions)
     ALLOCATE (elements(I)%ions(nions))
     ! Loop over all ions of given chem.element
     DO J = lowerion, upperion
        elements(I)%ions(J)%ion_stage = current_ion
        current_ion = lowerion + 1
     END DO
  END DO

  ! Only for testing
  PRINT*, 'testing'
  DO I = 1, n_elements
     element_index = I
     Z = elements(I)%atom_number
     lowerion = elements(I)%ions(1)%ion_stage
     upperion = elements(I)%ions(nions)%ion_stage
     PRINT*, element_index, Z, lowerion, upperion
  END DO

END SUBROUTINE read_composition


