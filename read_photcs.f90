SUBROUTINE read_photcs(inputdata, indexe, max_levels, phot_file)
USE types
USE constants
IMPLICIT NONE

! input variables
CHARACTER(LEN=20)                       :: phot_file
INTEGER                                 :: indexe, max_levels
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
INTEGER                                 :: l_index
DOUBLE PRECISION, PARAMETER             :: rydberg = 13.5979996 !(eV)
! reading files
INTEGER                                 :: ios
INTEGER                                 :: new_ion, old_ion
LOGICAL                                 :: save_cs
INTEGER                                 :: atom_number
INTEGER                                 :: electron_number

atom_number = elements(indexe)%atom_number

OPEN(13, FILE=phot_file)

SELECT CASE(inputdata)

! _____________________________________________________________________
! data from the Opacity project
! we expect datafile in this format
! *
! level index, atomic number, ion index, *, *, energy (Ryd), number of points
! frequency/treshold frequency, cross section () ! podívat se, v jakých je to jednotkách
CASE(2)
  ! a user can set maximal number of cross sectionw which will be read
  ! for the given ion, this is the reason why variables new_ion and old_ion
  ! are defined: we have to know, how many cross sections are read for the
  ! given ion
  old_ion = 0
  ! do 01: reading loop
  DO
   ! the first line contains information about the given excitation state
   READ(13, '(A)', IOSTAT=ios) line 
   ! write(*,*)  'read_photcs: ', line
   IF(line(1:1) == '*') CYCLE
   IF(ios /= 0) EXIT
   READ(line, *) l_index, indexZ, electron_number, junk, junk, energy, nofPoints
   IF(energy > 0.D0) THEN
    ! write(*,*) 'read_photcs: nofPoints = ', nofPoints
    DO I = 1, nofPoints
    END DO
   END IF
   ! we have to know, if the ion is different from the previous one
   new_ion = electron_number
   indexI = atom_number - electron_number + 1
   ! if it is different, we have to set the new variables
   IF(new_ion /= old_ion) THEN
    n_read = 0
    save_cs = .TRUE.
    old_ion = new_ion
    !print*, 'indexclev = ', indexclev
   END IF
   !print*, 'indexclev = ', indexclev, ' lowering_index = ', lowering_index
   n_read = n_read + 1
   IF(n_read > max_levels .AND. max_levels > 0) save_cs = .FALSE.
   IF(save_cs .EQV. .TRUE.) THEN
    CALL find_photion_elindex(indexe, indexI, l_index, indexclev)
    IF(indexclev < 0) THEN
     DO I = 1, nofPoints
      READ(13, '(A)', IOSTAT=ios) line 
     END DO
     CYCLE
    END IF
    ! write(*,*) 'read_photcs: indexclev = ', indexclev
    ALLOCATE(elements(indexe)%ions(indexI)%levels(indexclev)%photcros(2,nofPoints))
    IF(nofPoints /= 0) n_photcrossect = n_photcrossect + 1
    !write(*,*) 'read_photcs: n_photcrossect = ', n_photcrossect, 'nofPoints = ', nofPoints
    ! we can compute a frequency treshold from these data
    freqt = (MINVAL(elements(indexe)%ions(indexI + 1)%levels(:)%exci_energy) - &
     elements(indexe)%ions(indexI)%levels(indexclev)%exci_energy) / const_h
    elements(indexe)%ions(indexI)%levels(indexclev)%phfreq = freqt
   END IF
   !print*, 'nofPoints = ', nofPoints
   ! write(*,*) 'read_photcs: freqt = ', freqt
   !elements(indexe)%ions(indexI)%levels(indexclev)%phfreq = freqt
   ! now we will read the given data for the photoionization cross section
   ! write(*,*) 'read_photcs: nofPoints = ', nofPoints
   DO I = 1, nofPoints
    READ(13, '(A)', IOSTAT=ios) line
    if(ios /= 0) EXIT
    ! write(*,*) 'read_photcs: ', line
    READ(line,*) freq, cross
    !print*, I, freq, cross
    !print*, freq(I), cross(I)
    IF(save_cs .EQV. .TRUE.) THEN
     !freq = freq * freqt
     ! elements(indexe)%ions(indexI)%levels(indexclev)%photcros(1,I) = freq * freqt
     elements(indexe)%ions(indexI)%levels(indexclev)%photcros(1,I) = freq * Rydberg * e_v / const_h
     ! write(*,*) 'read_photcs: photfreq = ', freq * Rydberg * e_v / h
     elements(indexe)%ions(indexI)%levels(indexclev)%photcros(2,I) = cross * 1.D-18
     ! write(32,*) indexe, indexI, indexclev, freq * Rydberg * e_v / h, cross * 1.D-15
     ! write(32,*) indexe, indexI, indexclev, freq, cross
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
