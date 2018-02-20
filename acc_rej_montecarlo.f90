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
! INTEGER                       :: ios
DOUBLE PRECISION              :: sinseed, cosseed
! the linear interpolation parameters
DOUBLE PRECISION              :: a_linint, b_linint
DOUBLE PRECISION              :: x1, x2, fx1, fx2
! division on several intervals
DOUBLE PRECISION                :: summ, summ2
DOUBLE PRECISION                :: tot_int_val, tot_int_val2
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
! write(*,*) 'acc_rej_montecarlo: NR = ', NR
freq_min = incomingflux(2,1)
freq_max = incomingflux(NR, 1)
!print*, 'incomingflux(:,2): ', incomingflux(:,2)
flux_max = MAXVAL(incomingflux(:,2))
!print*, 'MAXVAL? ', junk
!print*, 'initial calculations for frequency generation: ', freq_min, freq_max, flux_max
! speed up a procedure: divide the whole file into several frequency intervals
! dividing the whole input file onto N intervals
IF(NR > 1000) THEN
 Nints = 100
ELSE IF(NR <=1000) THEN
 Nints = 50
END IF
ALLOCATE(int_val(Nints))
! the edges are now (for computing an index one has to add 1)
! I_1 = [1, N]
! I_2 = [N, 2N]
! .
! .
! .
! I_n = [n(N - 1), N (n + 1)]
! .
! .
! .
! I_N = [N(N - 1), NR]
Npoints = NR / Nints
! write(*,*) 'acc_rej_montecarlo: Npoints = ', Npoints, ' NR = ', NR
DO I = 1, Nints
 low_bound = 1 + Npoints * (I - 1)
 upp_bound = 1 + Npoints * I
 IF(I == Nints) upp_bound = NR
 summ = 0.D0
 DO J = low_bound, upp_bound - 1
  summ = summ + 5.D-1 * (incomingflux(J + 1, 2) + incomingflux(J, 2)) * &
   (incomingflux(J + 1, 1) - incomingflux(J, 1)) 
 END DO
 int_val(I) = summ
 ! write(*,*) 'acc_rej_montecarlo: low = ', low_bound, ' up = ', upp_bound, ' int_val = ', int_val(I)
END DO

! the total area under the graph calculation
summ = 0.D0
DO I = 1, Nints
 summ = summ + int_val(I)  
END DO
tot_int_val = summ
! summ2 = 0.D0
! DO I = 1, NR - 1
!   summ2 = summ2 + 5.D-1 * (incomingflux(I + 1, 2) + incomingflux(I, 2)) * &
!    (incomingflux(I + 1, 1) - incomingflux(I, 1)) 
! END DO
! tot_int_val2 = summ2
! write(*,*) 'acc_rej_montecarlo: tot_int_val = ', tot_int_val, ' tot_int_val2 = ', tot_int_val2
! generation of the photonic frequency with the given distribution
n_created_packs = 0
DO I = 1, Nints
 low_bound = 1 + Npoints * (I - 1)
 upp_bound = 1 + Npoints * I
 IF(I == Nints) upp_bound = NR
 loc_npacks = int_val(I) / tot_int_val * DBLE(n_packs)
 !write(*,*) 'acc_rej_montecarlo: loc_npacks = ', loc_npacks
 IF(loc_npacks == 0) CYCLE
 ! we generate now loc_npacks on the given interval
 DO K=n_created_packs + 1, n_created_packs + loc_npacks + 1
  FOUND = .FALSE.
  loc_freq_min = incomingflux(low_bound, 1)
  loc_freq_max = incomingflux(upp_bound, 1)
  loc_flux_max = MAXVAL(incomingflux(low_bound:upp_bound, 2))
  DO WHILE ( FOUND .EQV. .FALSE.)
   CALL random_number(sinseed)
   CALL random_number(cosseed)
   !sinseed=INT(100*abs(sin(REAL(n_pack))))
   !cosseed=INT(100*abs(cos(seed)))
   ran_freq = loc_freq_min + (loc_freq_max - loc_freq_min) * sinseed
   ran_flux = cosseed * flux_max
   !print*, 'random frequency and flux values: ', ran_freq, ran_flux
   ! now we have to say if this points are located under the flux curve
   ! how we do it? 
   ! 1.) find if the random frequency is on some point in the frequency file
   ! then the bound is readable directly from the file
   ! 2.) if the random frequency is between two points we have to interpolate
   ! between these two flux points and then calculate the given bound point
   ! 1.)
    !print*, 'and we will calculate the linear interpolation...'
   ! i find two frequency points between we will interpolate
   DO J=low_bound, upp_bound
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
   ! print*, 'lw_index = ', lw_index, 'ran_freq = ', ran_freq, 'ran_flux', ran_flux, ' bound_flux = ', bound_flux
   IF (ran_flux < bound_flux) THEN
    freq(K) = ran_freq
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
 n_created_packs = n_created_packs + loc_npacks
END DO
write(*,*) 'acc_rej_montecarlo: n_created_packs = ', n_created_packs, ' n_pack = ', n_packs


DO K=n_created_packs + 1, n_packs
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
  IF (FOUND .EQV. .TRUE.) THEN
  ! print*, 'we found the photon :-)'
  ELSE
  ! print*, 'we did not find it ... trying again...'
  END IF
  ! if ((MODULO(K,10000) .EQ. 0) .AND. (FOUND .EQV. .TRUE.)) print*, 'Generating the frequency packet ', K, ' ...'
 END DO
END DO
OPEN(37, FILE='freq_dist_test.dat')
 DO I = 1, n_packs
  WRITE(37, *) freq(I)
 END DO
CLOSE(37)
! STOP 'acc_rej_montecarlo, testing'
 ! after the last photon is calculated we erase the field incomingflux'
! IF(n_packet .EQ. n_packs) DEALLOCATE(incomingflux)
END SUBROUTINE acc_rej_montecarlo
