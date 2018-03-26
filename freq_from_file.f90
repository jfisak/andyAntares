! this subroutine generates a random frequency of photon
! using distribution from an emergent flux computed by
! another model of atmosphere
! this procedure is called if and only if incomingflux = 1
SUBROUTINE freq_from_file(n_packs,freq)

  USE types

  IMPLICIT NONE 

  INTEGER                       :: NR,n_packet,n_packs
  INTEGER                       :: I,J, lw_index, lg_index
!  LOGICAL                       :: found, linint
  INTEGER, PARAMETER            :: maxrows = 6000000
  DOUBLE PRECISION              :: seed
  DOUBLE PRECISION, DIMENSION(n_packs) :: freq
!  DOUBLE PRECISION              :: freq, freq_min, freq_max, flux_max
!  DOUBLE PRECISION              :: ran_freq, ran_flux, bound_flux
  DOUBLE PRECISION              :: junk
  INTEGER                       :: ios
! (2) PoWR testing model
DOUBLE PRECISION                        :: logwv
CHARACTER(100)                  :: fluxfile
!  DOUBLE PRECISION              :: sinseed, cosseed
!  ! the linear interpolation parameters
!  DOUBLE PRECISION              :: a_linint, b_linint
!  DOUBLE PRECISION              :: x1, x2, fx1, fx2
 
! in the first photon computes the flux field which depends on frequency
! flux(NUMBER OF ROW, INDEX) INDEX = 1 ... FREQUENCY, INDEX = 2 ... FLUX
! at first we will use the given input flux
SELECT CASE (inputflux)
 !______________________________________________________________________________
 ! (1) TLUSTY model
 ! in this case we read input file from the TLUSTY model
 ! number   frequency
CASE (1)
 IF(ALLOCATED(incomingflux) .EQV. .FALSE.) THEN
   !print*, 'initialization of reading input flux'
   OPEN(11,status='old',FILE='emflux.dat')
   ! it is necessary to compute number of rows of the file
   NR = 0
  DO 
    READ(11,*,IOSTAT=ios) junk, junk
   IF (ios /= 0) EXIT
   NR = NR + 1
  END DO
  write(*,*) 'freq_from_file: nr = ', NR
  REWIND(11)
  ! now we can allocate the field flux (frequency, flux))
  write(*,*) 'allocation of the field incomingflux(', NR, ', 2)'
  ALLOCATE(incomingflux(NR,2))
  ! and read from given file
  PRINT*, 'reading the flux from input file...'
  DO I=1,NR
   READ(11,*) incomingflux(I,1), incomingflux(I,2)
  END DO
  CLOSE(11)
  do I=1,NR
  ! print*, 'incomingflux: ', incomingflux(I,1), ', ', incomingflux(I,2)
  end do
 END IF
 !______________________________________________________________________________
 ! (2) PoWR model
 ! log(wavelength)      Intensity
 ! log(Angstroms)       erg / cm^2 / s / Hz
 CASE(2)
  CALL GET_ENVIRONMENT_VARIABLE("FLUXFILE", fluxfile)
  IF(TRIM(fluxfile) == '') STOP 'freq_from_file: file was not found'
   OPEN(11, status='old', FILE=TRIM(fluxFile))
   IF(.NOT. ALLOCATED(incomingflux)) THEN
   NR = 0
   DO
    READ(11,*,IOSTAT=ios) junk, junk
    IF (ios /= 0) EXIT
    NR = NR + 1
   END DO
   ALLOCATE(incomingflux(NR, 2))
   REWIND(11)
   PRINT*, 'reading the flux from input file...'
   write(*,*) 'freq_from_file: NR = ', NR
   DO I=1,NR
    ! wl in log(Angstroms), flux in erg / cm^2 / s / Hz
    READ(11,*) logwv, incomingflux(I,2)
    ! write(*,*) 'freq_from_file: I = ', I, ' logwv = ', logwv, ' I = ', incomingflux(I, 2)
    incomingflux(I, 1) = light_speed / (1.D-10 * exp(logwv))
   END DO
   CLOSE(11)
  END IF
 CASE DEFAULT
  print*, 'the choice of variable inputflux = ', inputflux, 'is not known...'
  STOP 'ENDING PROGRAM NOW...'
END SELECT

CALL acc_rej_montecarlo(n_packs,freq)
  if ((MODULO(n_packet,10000) .EQ. 0)) print*, 'Generating the frequency packet ', n_packet, ' ...'
 IF(n_packet .EQ. n_packs) DEALLOCATE(incomingflux)
END SUBROUTINE freq_from_file
