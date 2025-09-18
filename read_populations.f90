! MUST BE COMPLETED!!!
! this sbr reads populations listed in a selected file
SUBROUTINE read_populations()

USE types
USE constants

IMPLICIT NONE

CHARACTER(filename_length)                             :: populationfile

populationfile=TRIM(inputpopFile)

write(*,*) 'read_populations: populationfile = ', populationfile

! allocation of all population arrays
OPEN(UNIT=13, FILE=populationfile)
 DO
  ! READ(8,'(A)') line
  ! IF(line(1:1) .EQ. '*') CYCLE
   
 END DO
CLOSE(13)

END SUBROUTINE read_populations
