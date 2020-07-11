! this subroutine generates a random frequency of photon
! using distribution from an emergent flux computed by
! another model of atmosphere
! this procedure is called if and only if incomingflux = 1
SUBROUTINE freq_from_file(n_packs,freq)

  USE types

  IMPLICIT NONE 

  INTEGER                       :: NR,n_packet,n_packs
  INTEGER                       :: I
!  LOGICAL                       :: found, linint
  INTEGER, PARAMETER            :: maxrows = 6000000
  DOUBLE PRECISION, DIMENSION(n_packs) :: freq
!  DOUBLE PRECISION              :: freq, freq_min, freq_max, flux_max
!  DOUBLE PRECISION              :: ran_freq, ran_flux, bound_flux
  DOUBLE PRECISION              :: junk
  INTEGER                       :: ios
! (2) PoWR testing model
DOUBLE PRECISION                        :: logwv, Iflux
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
  ! write(99,*) 'incomingflux: ', incomingflux(I,1), ', ', incomingflux(I,2)
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
    READ(11,*) logwv, Iflux
    ! CHECK ONCE MORE !!!!!!!!!!!!!!!!!
    incomingflux(NR - I + 1, 1) = 1.D6 * light_speed * exp(- logwv)
    ! incomingflux(I, 2) = light_speed * Iflux / (incomingflux(I, 1) ** 2.0)
    incomingflux(I, 2) = Iflux
    ! write(*,*) 'freq_from_file: I = ', I, ' logwv = ', logwv, ' freq = ', incomingflux(I, 1), ' flux = ', incomingflux(I, 2)
   END DO
   CLOSE(11)
  END IF
 CASE DEFAULT
  write(99,*) 'the choice of variable inputflux = ', inputflux, 'is not known...'
  STOP 'ENDING PROGRAM NOW...'
END SELECT

CALL acc_rej_montecarlo(n_packs,freq)
  if ((MODULO(n_packet,10000) .EQ. 0)) write(99,*) 'Generating the frequency packet ', n_packet, ' ...'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we have the field incomingflux allocated and defined values...now we can
! generate the photonic frequencies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we generate a new photonic frequency
! we will consider that frequencies are in the right order
! NR = INT(incomingflux(1,1))
! freq_min = incomingflux(2,1)
! freq_max = incomingflux(NR + 1, 1)
! !print*, 'incomingflux(:,2): ', incomingflux(:,2)
! flux_max = MAXVAL(incomingflux(:,2))
! !print*, 'MAXVAL? ', junk
! !print*, 'initial calculations for frequency generation: ', freq_min, freq_max, flux_max
! FOUND = .FALSE.
! DO WHILE ( FOUND .EQV. .FALSE.)
!  CALL random_number(sinseed)
!  CALL random_number(cosseed)
!  !sinseed=INT(100*abs(sin(REAL(n_pack))))
!  !cosseed=INT(100*abs(cos(seed)))
!  ran_freq = freq_min + (freq_max-freq_min)*sinseed
!  ran_flux = cosseed*flux_max
!  !print*, 'random frequency and flux values: ', ran_freq, ran_flux
!  ! now we have to say if this points are located under the flux curve
!  ! how we do it? 
!  ! 1.) find if the random frequency is on some point in the frequency file
!  ! then the bound is readable directly from the file
!  ! 2.) if the random frequency is between two points we have to interpolate
!  ! between these two flux points and then calculate the given bound point
!  ! 1.)
!  linint = .true.
!  DO I=2,NR+1
!   IF ( ran_freq .EQ. incomingflux(I,1)) THEN
!    bound_flux = incomingflux(I,2)
!    linint = .false.
!    !print*, 'linit is now in false...'
!   END IF
!  END DO
!  ! 2.)
!  IF ( linint .EQV. .true.) THEN
!    !print*, 'and we will calculate the linear interpolation...'
!   ! i find two frequency points between we will interpolate
!   DO J=2,NR+1
!    IF (ran_freq > incomingflux(J,1)) THEN
!     !print*, 'loop ', J, 'was cycled...'
!     CYCLE
!    ELSE
!     lw_index = J - 1
!     lg_index = J
!     !print*, 'exiting cycle in loop ', J
!     EXIT
!    END IF
!   END DO
!   ! ii linear interpolation
!    x1 = incomingflux(lw_index, 1)
!    x2 = incomingflux(lg_index, 1)
!    fx1 = incomingflux(lw_index, 2)
!    fx2 = incomingflux(lg_index, 2)
!    a_linint = (fx2 - fx1)/(x2 - x1)
!    b_linint = (fx1*x2 - fx2*x1)/(x2 - x1) 
!    ! now we can find bound flux easily
!    bound_flux = a_linint * ran_freq + b_linint
!   ! iii did we find the right frequency?
!   !print*, 'lw_index = ', lw_index, 'ran_freq = ', ran_freq, 'ran_flux', ran_flux, ' bound_flux = ', bound_flux
!   IF (ran_flux < bound_flux) THEN
!    freq = ran_freq
!    FOUND = .TRUE.
!   END IF
!  END IF
!  IF (FOUND .EQV. .TRUE.) THEN
!  ! print*, 'we found the photon :-)'
!  ELSE
!  ! print*, 'we did not find it ... trying again...'
!  END IF
! END DO
 ! after the last photon is calculated we erase the field incomingflux'
 IF(n_packet .EQ. n_packs) DEALLOCATE(incomingflux)
END SUBROUTINE freq_from_file
