! reads populations from a file
!
! INPUT: fileFound(CHAR(filename_length)) -- name of file
! OUTPUT: NONE
!
SUBROUTINE find_populations(fileFound)

USE types
USE constants

IMPLICIT NONE

LOGICAL                                 :: fileFound
CHARACTER(LEN=filename_length)                       :: populFile, inputfile
INTEGER                                 :: num_rows, npop_rows
! reading from file
CHARACTER(LEN=100)                      :: line
INTEGER                                 :: reading_populations
CHARACTER(LEN=filename_length)                       :: popInpFile
INTEGER                                 :: indexe, indexi, part_Z
INTEGER                                 :: dataType
! type for a better data saving
TYPE inp_Ion
 INTEGER                                :: indexe, indexi, inpType
 CHARACTER(LEN=filename_length)                      :: fileName
END TYPE inp_Ion
TYPE(inp_Ion), ALLOCATABLE               :: inpIon(:)
INTEGER                                 :: ind_I, gridcell
DOUBLE PRECISION                        :: logntot, loggi
! saving populations
LOGICAL                                 :: fileExists


CALL GET_ENVIRONMENT_VARIABLE("POPULATIONS", populFile)

! write(*,*) 'find_populations: populFile = ', populFile
IF(populFile == '') THEN
 ! no input data were found, the data will be calculated by the code
 fileFound = .FALSE.
 RETURN
ELSE
 ! the code will read the data which will be used for the first iteration
 fileFound = .TRUE.
END IF

write(inputfile,"(A, A)") trim(outputfolder), trim(populFile)

! user can set the variable POPULATIONS but the file still could not exist
inquire( file=inputfile, exist=fileFound )

IF(.NOT. fileFound) RETURN

OPEN(51,FILE=inputfile)
 num_rows = 0
 DO 
  READ(51,'(A)', IOSTAT = reading_populations) line
  IF(reading_populations /= 0) EXIT
  IF(line(1:1) .EQ. '*') CYCLE
  num_rows = num_rows + 1
 END DO
 REWIND(51)
 ALLOCATE(inpIon(num_rows))
 ind_I = 0
 line = ' '
 DO
  READ(51,'(A)', IOSTAT = reading_populations) line
  IF(reading_populations /= 0) EXIT
  ! write(*,*) 'find_populations: #1 line = ', line
  IF(line(1:1) .EQ. '*') CYCLE
  READ(line,*) part_Z, indexi, dataType, popInpFile
  ! write(*,*) 'find_populations: #2 popInpFile = ', popInpFile, LEN(popInpFile)
  ! write(*,*) 'find_populations: #3 Z = ', Z, ' indexi = ', indexi,&
  !  'dataType = ', dataType 
  ! write(*,*) 'find_populations: #4 popInpFile = ', popInpFile, LEN(popInpFile)
  ! IF(Z == 1 .AND. indexi == 1) popInpFile = '/home/jakub/Documents/PhD/3dwind/PoWRtestCase/H_I.dat'
  ! IF(Z == 1 .AND. indexi == 2) popInpFile = '/home/jakub/Documents/PhD/3dwind/PoWRtestCase/H_II.dat'
  ! IF(Z == 2 .AND. indexi == 1) popInpFile = '/home/jakub/Documents/PhD/3dwind/PoWRtestCase/He_I.dat'
  ! IF(Z == 2 .AND. indexi == 2) popInpFile = '/home/jakub/Documents/PhD/3dwind/PoWRtestCase/He_II.dat'
  ! IF(Z == 2 .AND. indexi == 3) popInpFile = '/home/jakub/Documents/PhD/3dwind/PoWRtestCase/He_III.dat'
  write(*,*) 'find_populations: part_Z = ', part_Z, ' indexi = ', indexi,&
   'dataType = ', dataType, ' popInpFile = ', TRIM(popInpFile)
  ind_I = ind_I + 1
  inpIon(ind_I)%indexe = part_Z
  inpIon(ind_I)%indexi = indexi
  inpIon(ind_I)%fileName = popInpFile
  inpIon(ind_I)%inpType = dataType
 END DO
CLOSE(51)

! reading the populations from the given files
DO ind_I = 1, num_rows
 indexe = inpIon(ind_I)%indexe
 indexi = inpIon(ind_I)%indexi
 popInpFile = inpIon(ind_I)%fileName
 dataType = inpIon(ind_I)%inpType
 ! write(*,*) 'find_populations: indexe = ', indexe, ' indexi = ', indexi
 SELECT CASE(dataType)
  ! PoWR test case
  CASE(3)
   inquire( file=popInpFile, exist=fileExists )
   IF(.NOT. fileExists) CYCLE
   OPEN(52, FILE=popInpFile)
   npop_rows = 0
   DO 
    READ(52,'(A)', IOSTAT = reading_populations) line
    IF(reading_populations /= 0) EXIT
    IF(line(1:1) .EQ. '*') CYCLE
    npop_rows = npop_rows + 1
   END DO
   REWIND(52)
   DO gridcell = 1, npop_rows
    READ(52,'(A)',IOSTAT=reading_populations) line
    ! write(*,*) 'find_populations: ', line
    IF(reading_populations /= 0) EXIT
    READ(line,*) logntot, loggi
    model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop = 10**(loggi+logntot)
    model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = 10**logntot
    ! write(*,*) 'find_populations: gridcell = ', gridcell, ' indexe = ', indexe, ' indexi = ', indexi
    ! write(*,'(A,E12.3,A,E12.3)') 'find_populations: gl_pop = ', 10**(loggi+logntot), ' tot_pop = ', 10**logntot
   END DO
   CLOSE(52)
  CASE DEFAULT
   write(99,*) 'find_populations: wrong populations type'
   write(99,*) ' in the file ', populFile
   STOP
 END SELECT
END DO


END SUBROUTINE find_populations
