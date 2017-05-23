SUBROUTINE read_photcs(inputdata, element, max_levels, filename)
USE types
IMPLICIT NONE

! input variables
CHARACTER(LEN=20)                       :: filename
INTEGER                                 :: element, ion, max_levels
INTEGER                                 :: inputdata
CHARACTER(LEN=20)                       :: junk
! loop index
INTEGER                                 :: I
! informations about element
INTEGER                                 :: elementI, elementZ, indexclev
DOUBLE PRECISION                        :: energy
! treshold frequency
INTEGER                                 :: nofPoints
DOUBLE PRECISION                        :: freqt
! data about cross section
DOUBLE PRECISION, ALLOCATABLE           :: freq(:), cross(:)
INTEGER                                 :: n_read
DOUBLE PRECISION, PARAMETER             :: rydberg = 13.5979996 !(eV)
! reading files
INTEGER                                 :: ios


SELECT CASE(inputdata)

! _____________________________________________________________________
! data from the Opacity project
! we expect datafile in this format
! *
CASE(2)
 OPEN(9, FILE=filename)
  n_read = 0
  ! do 01: reading loop
  DO
  ! the first line contains information about the given excitation state
   READ(9,IOSTAT=ios) indexclev, elementZ, elementI, junk, junk, energy, nofPoints
   IF(ios /= 0) EXIT
   ALLOCATE(freq(nofPoints), cross(nofPoints), &
        elements(element)%ions(elementI)%levels(indexclev)%photcros(1,nofPoints))
   ! we can compute a frequency treshold from these data
   freqt = abs(energy * Rydberg) / h
   ! now we will read the given data for the photoionization cross section
   DO I = 1, nofPoints
    READ(9) freq(I), cross(I)
    freq(I) = freq(I) * freqt
    elements(element)%ions(elementI)%levels(indexclev)%photcros(1,I) = freq(I)
    elements(element)%ions(elementI)%levels(indexclev)%photcros(2,I) = cross(I)
   END DO
   ! it would be possible to save these data in the form of fit coefficients
   ! or to save data itself (I don't know if only fit is good enough)
   ! we need to deallocate these fields before reading another data
   DEALLOCATE(freq, cross)
   n_read = n_read + 1
   IF(n_read == max_levels) EXIT
  ! end do 01: reading loop
  END DO
CASE DEFAULT
 STOP 'this input for photoionization cross sections is not known...'
END SELECT

END SUBROUTINE read_photcs
