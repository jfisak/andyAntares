! this procedure will choose a frequency using
! acception/rejection method Monte Carlo
SUBROUTINE acc_rej_montecarlo(n_packs,freq)

  USE types

  IMPLICIT NONE 

  INTEGER                       :: NR,n_packet,n_packs
  INTEGER                       :: I,J, K, lw_index, lg_index
  LOGICAL                       :: found, linint
  INTEGER, PARAMETER            :: maxrows = 6000000
  DOUBLE PRECISION              :: seed
  DOUBLE PRECISION, DIMENSION(n_packs) :: freq
  DOUBLE PRECISION              :: freq_min, freq_max, flux_max
  DOUBLE PRECISION              :: ran_freq, ran_flux, bound_flux
  DOUBLE PRECISION              :: junk
!  INTEGER                       :: ios
  DOUBLE PRECISION              :: sinseed, cosseed
  ! the linear interpolation parameters
  DOUBLE PRECISION              :: a_linint, b_linint
  DOUBLE PRECISION              :: x1, x2, fx1, fx2
 
! in the first photon computes the flux field which depends on frequency
! flux(NUMBER OF ROW, INDEX) INDEX = 1 ... FREQUENCY, INDEX = 2 ... FLUX
! at first we will use the given input flux
!SELECT CASE (inputflux)
! ! in this case we read input file from Jiri Kubat's model in the form
! ! row number   frequency       wavelength      flux    log flux
!CASE (1)
! IF(n_packet .EQ. 1) THEN
!   !print*, 'initialization of reading input flux'
!   OPEN(11,status='old',FILE='emflux.dat')
!   ! it is necessary to compute number of rows of the file
!   NR = 0
!  DO I=1,maxrows
!    READ(11,*,IOSTAT=ios) junk, junk, junk, junk, junk
!   IF (ios /= 0) EXIT
!   IF (I == maxrows) THEN
!    print*, 'Error: Maximum number of records exceeded...'
!    print*, 'Exiting program now...'
!    STOP
!   END IF
!   NR = NR + 1
!  END DO
!  REWIND(11)
!  ! now we can allocate the field flux (frequency, flux))
!  IF (DEBUG .EQ. 1) print*, 'allocation of the field incomingflux(', NR, ', 2)'
!  ALLOCATE(incomingflux(NR+1,2))
!  ! and read from given file
!  PRINT*, 'reading the flux from input file...'
!  incomingflux(1,1) = NR
!  incomingflux(1,2) = 0
!  DO I=2,NR+1
!   READ(11,*) junk, incomingflux(I,1), junk, incomingflux(I,2), junk
!  END DO
!  CLOSE(11)
!  do I=1,NR+1
!  ! print*, 'incomingflux: ', incomingflux(I,1), ', ', incomingflux(I,2)
!  end do
! END IF
! CASE DEFAULT
!  print*, 'the choice of variable inputflux = ', inputflux, 'is not known...'
!  STOP 'ENDING PROGRAM NOW...'
!END SELECT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we have the field incomingflux allocated and defined values...now we can
! generate the photonic frequencies
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we generate a new photonic frequency
! we will consider that frequencies are in the right order
 NR = INT(incomingflux(1,1))
 freq_min = incomingflux(2,1)
 freq_max = incomingflux(NR + 1, 1)
 !print*, 'incomingflux(:,2): ', incomingflux(:,2)
 flux_max = MAXVAL(incomingflux(:,2))
 !print*, 'MAXVAL? ', junk
 !print*, 'initial calculations for frequency generation: ', freq_min, freq_max, flux_max
DO K=1,n_packs
 FOUND = .FALSE.
 DO WHILE ( FOUND .EQV. .FALSE.)
  CALL random_number(sinseed)
  CALL random_number(cosseed)
  !sinseed=INT(100*abs(sin(REAL(n_pack))))
  !cosseed=INT(100*abs(cos(seed)))
  ran_freq = freq_min + (freq_max-freq_min)*sinseed
  ran_flux = cosseed*flux_max
  !print*, 'random frequency and flux values: ', ran_freq, ran_flux
  ! now we have to say if this points are located under the flux curve
  ! how we do it? 
  ! 1.) find if the random frequency is on some point in the frequency file
  ! then the bound is readable directly from the file
  ! 2.) if the random frequency is between two points we have to interpolate
  ! between these two flux points and then calculate the given bound point
  ! 1.)
  linint = .true.
  DO I=2,NR+1
   IF ( ran_freq .EQ. incomingflux(I,1)) THEN
    bound_flux = incomingflux(I,2)
    linint = .false.
    !print*, 'linit is now in false...'
   END IF
  END DO
  ! 2.)
  IF ( linint .EQV. .true.) THEN
    !print*, 'and we will calculate the linear interpolation...'
   ! i find two frequency points between we will interpolate
   DO J=2,NR+1
    IF (ran_freq > incomingflux(J,1)) THEN
     !print*, 'loop ', J, 'was cycled...'
     CYCLE
    ELSE
     lw_index = J - 1
     lg_index = J
     !print*, 'exiting cycle in loop ', J
     EXIT
    END IF
   END DO
   ! ii linear interpolation
    x1 = incomingflux(lw_index, 1)
    x2 = incomingflux(lg_index, 1)
    fx1 = incomingflux(lw_index, 2)
    fx2 = incomingflux(lg_index, 2)
    a_linint = (fx2 - fx1)/(x2 - x1)
    b_linint = (fx1*x2 - fx2*x1)/(x2 - x1) 
    ! now we can find bound flux easily
    bound_flux = a_linint * ran_freq + b_linint
   ! iii did we find the right frequency?
   !print*, 'lw_index = ', lw_index, 'ran_freq = ', ran_freq, 'ran_flux', ran_flux, ' bound_flux = ', bound_flux
   IF (ran_flux < bound_flux) THEN
    freq(K) = ran_freq
    FOUND = .TRUE.
   END IF
  END IF
  IF (FOUND .EQV. .TRUE.) THEN
  ! print*, 'we found the photon :-)'
  ELSE
  ! print*, 'we did not find it ... trying again...'
  END IF
  if ((MODULO(K,10000) .EQ. 0) .AND. (FOUND .EQV. .TRUE.)) print*, 'Generating the frequency packet ', K, ' ...'
 END DO
END DO
 ! after the last photon is calculated we erase the field incomingflux'
! IF(n_packet .EQ. n_packs) DEALLOCATE(incomingflux)
END SUBROUTINE acc_rej_montecarlo
