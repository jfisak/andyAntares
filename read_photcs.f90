SUBROUTINE read_photcs(inputdata, indexe, max_levels, phot_file)
USE types
IMPLICIT NONE

! input variables
CHARACTER(LEN=20)                       :: phot_file
INTEGER                                 :: indexe, ion, max_levels
INTEGER                                 :: inputdata
CHARACTER(LEN=20)                       :: junk
! loop index
CHARACTER(LEN=200)                      :: line
INTEGER                                 :: I
! informations about element
INTEGER                                 :: indexi, indexZ, indexclev
DOUBLE PRECISION                        :: energy
! treshold frequency
INTEGER                                 :: nofPoints
DOUBLE PRECISION                        :: freqt
! data about cross section
DOUBLE PRECISION                        :: freq, cross
INTEGER                                 :: n_read
DOUBLE PRECISION, PARAMETER             :: rydberg = 13.5979996 !(eV)
! energy calculation
DOUBLE PRECISION                        :: ionstage, energyFZ
! reading files
INTEGER                                 :: ios
INTEGER                                 :: new_ion, old_ion
INTEGER                                 :: lowering_index
LOGICAL                                 :: save_cs

OPEN(13, FILE=phot_file)

SELECT CASE(inputdata)

! _____________________________________________________________________
! data from the Opacity project
! we expect datafile in this format
! *
CASE(2)
  old_ion = 0
  ! do 01: reading loop
  DO
   ! the first line contains information about the given excitation state
   READ(13, '(A)', IOSTAT=ios) line
   !print*, 'read_photcs: ', line
   IF(ios /= 0) EXIT
   READ(line, *) indexclev, indexZ, indexI, junk, junk, energy, nofPoints
   !print*, 'read_photcs: indexclev, indexZ, nofPoints', indexclev, indexZ, nofPoints
   ! calculation of total energy
!   ionstage = 0.D0
!   DO I = 1, indexI
!    ionstage = ionstage + elements(index)%ions(indexI)%ion_potential
!   END DO
!   energyFZ = rydberg * energy * e_v + ionstage
!   print*, 'read_photcs: energyFZ = ', energyFZ / e_v / rydberg
   new_ion = indexI
   IF(new_ion /= old_ion) THEN
    n_read = 0
    save_cs = .TRUE.
    old_ion = new_ion
    lowering_index = indexclev
    !print*, 'indexclev = ', indexclev
   END IF
   indexclev = indexclev - lowering_index + 1
   !print*, 'indexclev = ', indexclev, ' lowering_index = ', lowering_index
   n_read = n_read + 1
   IF(n_read > max_levels .AND. max_levels > 0) save_cs = .FALSE.
   IF(save_cs .EQV. .TRUE.) THEN
    ALLOCATE(elements(indexe)%ions(indexI)%levels(indexclev)%photcros(2,nofPoints))
    n_photcrossect = n_photcrossect + 1
   END IF
   !print*, 'nofPoints = ', nofPoints
   ! we can compute a frequency treshold from these data
   freqt = abs(energy * Rydberg * e_v) / h
   !elements(indexe)%ions(indexI)%levels(indexclev)%phfreq = freqt
   ! now we will read the given data for the photoionization cross section
   DO I = 1, nofPoints
    READ(13, '(A)', IOSTAT=ios) line
    if(ios /= 0) EXIT
    !print*, 'read_photcs: ', line
    READ(line,*) freq, cross
    !print*, I, freq, cross
    !print*, freq(I), cross(I)
    IF(save_cs .EQV. .TRUE.) THEN
     !freq = freq * freqt
     elements(indexe)%ions(indexI)%levels(indexclev)%photcros(1,I) = freq * freqt
     elements(indexe)%ions(indexI)%levels(indexclev)%photcros(2,I) = cross * 1.D-15
     ! print*, index, indexI, freq(I), cross(I)
    END IF
   END DO
   ! it would be possible to save these data in the form of fit coefficients
   ! or to save data itself (I don't know if only fit is good enough)
   ! we need to deallocate these fields before reading another data
  ! end do 01: reading loop
  END DO
CASE DEFAULT
 STOP 'this input for photoionization cross sections is not known...'
END SELECT

CLOSE(13)
END SUBROUTINE read_photcs
