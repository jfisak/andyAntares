! this procedure will choose a frequency using
! acception/rejection method Monte Carlo
SUBROUTINE acc_rej_montecarlo(n_packs,freq)

USE types

IMPLICIT NONE 

INTEGER                       :: NR, n_packs
INTEGER                       :: I,J, K, lw_index, lg_index
LOGICAL                       :: found
INTEGER, PARAMETER            :: maxrows = 6000000
DOUBLE PRECISION, DIMENSION(n_packs) :: freq
DOUBLE PRECISION              :: freq_min, freq_max, flux_max
DOUBLE PRECISION              :: ran_freq, ran_flux, bound_flux
! INTEGER                       :: ios
DOUBLE PRECISION              :: sinseed, cosseed
! the linear interpolation parameters
DOUBLE PRECISION              :: a_linint, b_linint
DOUBLE PRECISION              :: x1, x2, fx1, fx2
! division on several intervals
DOUBLE PRECISION                :: summ
DOUBLE PRECISION                :: tot_int_val
INTEGER                         :: low_bound, upp_bound
INTEGER                         :: Nints, Npoints
DOUBLE PRECISION, ALLOCATABLE   :: int_val(:)
! local number of photons
INTEGER                         :: loc_npacks, n_created_packs
DOUBLE PRECISION                :: loc_flux_max, loc_freq_min, loc_freq_max
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! now we generate a new photonic frequency
! we will consider that frequencies are in the right order
NR = SIZE(incomingflux(:,1))
freq_min = incomingflux(1,1)
freq_max = incomingflux(NR, 1)
flux_max = MAXVAL(incomingflux(:,2))
! speed up a procedure: divide the whole file into several frequency intervals
! dividing the whole input file onto N intervals
IF(NR > 1000) THEN
 Nints = 100
ELSE IF(NR <=1000) THEN
 Nints = 50
END IF

DO K = 1, n_packs
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
  ! print*, 'and we will calculate the linear interpolation...'
  ! i find two frequency points between we will interpolate
  DO J=1, NR
   IF (ran_freq < incomingflux(J,1)) THEN
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
  IF (ran_flux < bound_flux) THEN
   freq(K) = ran_freq
   write(39,*) ran_freq
   FOUND = .TRUE.
  END IF
  IF (FOUND .EQV. .TRUE.) THEN
  ! print*, 'we found the photon :-)'
  ELSE
  ! print*, 'we did not find it ... trying again...'
  END IF
  ! if ((MODULO(K,10000) .EQ. 0) .AND. (FOUND .EQV. .TRUE.)) print*, 'Generating the frequency packet ', K, ' ...'
 END DO
END DO
! OPEN(37, FILE='freq_dist_test.dat')
!  DO I = 1, n_packs
!   WRITE(37, *) freq(I)
!  END DO
! CLOSE(37)
! STOP 'acc_rej_montecarlo, testing'
 ! after the last photon is calculated we erase the field incomingflux'
! IF(n_packet .EQ. n_packs) DEALLOCATE(incomingflux)
END SUBROUTINE acc_rej_montecarlo
