SUBROUTINE resonance_distance(pack_index, f_line, cell_dist, ldist, inCell)
USE types
USE counters
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index
DOUBLE PRECISION                :: cell_dist
DOUBLE PRECISION                :: f_line
! output variables
DOUBLE PRECISION                :: ldist
! beta law 
LOGICAL                         :: active
DOUBLE PRECISION, DIMENSION(3)  :: halfpos, lowbond, upbond, rbound
! DOUBLE PRECISION, PARAMETER     :: minint = 1.D-15
DOUBLE PRECISION                :: minint
DOUBLE PRECISION                :: chint
! testing
DOUBLE PRECISION                :: ldist_analyt
INTEGER                         :: n_pack_d, dummypackage
DOUBLE PRECISION                :: bfreq, halffreq, D
INTEGER                         :: OMP_GET_THREAD_NUM
DOUBLE PRECISION                :: ufreq, lfreq
LOGICAL                         :: TESTING = .TRUE.
LOGICAL                         :: inCell
LOGICAL                         :: calculate
INTEGER                         :: I
INTEGER                         :: cell_number
DOUBLE PRECISION                :: dist
DOUBLE PRECISION                :: f_res

dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
! write(*,*) 'resonance_distance: pack_index = ', pack_index, ' f_line = ', f_line,&
!  ' cell_dist = ', cell_dist

! write(*,*) '*****************************************************'
! 
! write(*,*) 'resonance_distance: fcmf / f_line = ', package(pack_index)%freq_cmf / f_line

! write(*,*) 'resonance_distance***************************************'
minint = 1.D0 / linelist(1)%freq
calculate = .TRUE.
! minint = 1.D-4
IF(package(pack_index)%freq_cmf < f_line) THEN
 calculate = .FALSE.
END IF
! Calculate distance the photon needs to travel to come to
! resonance with the next line. This assumes homologous
! expansion i.e. velocity is proportional to r. Projected
! gradient of the projected velocity in a direction of the
! photon propagation is more complicated in case of no
! homologous expansion
! homologous approximation
IF(TESTING .EQV. .TRUE.) THEN
 ldist_analyt = light_speed * ( R_inf / V_inf ) * &
  ( ( package(pack_index)%freq_cmf / package(pack_index)%freq_rf)-&
   f_line / package(pack_index)%freq_rf)
 package(dummypackage) = package(pack_index)
 inCell = .TRUE.
 ldist = ldist_analyt
 RETURN
 ! write(*,*) 'resonance_distance: ldist_analyt = ', ldist_analyt / R_inf
END IF
 ! boundary coordinate
! write(*,*) 'resonance_distance: cell_dist = ', cell_dist
CALL move_package(dummypackage, cell_dist)
bfreq = package(dummypackage)%freq_cmf
rbound = package(dummypackage)%pos
! frequency in the propagation cell boundary
! package(dummypackage)%pos = rbound
lowbond = package(pack_index)%pos
lfreq = package(pack_index)%freq_cmf
upbond = rbound
! cmf frequency on boundary
IF(lfreq < f_line) THEN
 calculate = .FALSE.
END IF
ufreq = bfreq
! we can calculate the line_dist
! lower and upper boundary
! write(*,*) 'resonance_distance: bfreq / f_line = ', bfreq / f_line
IF(bfreq < f_line .AND. calculate .EQV. .TRUE.) THEN
 active = .TRUE.
 inCell = .TRUE.
 ! write(*,*) 'resonance_distance: line res in prop cell'
 I = 0
 DO WHILE(active .EQV. .TRUE.)
  I = I + 1
  ! write(*,*) 'resonance_distance: I = ', I
  IF(I == 1000) THEN
   package(pack_index)%active = 0
   ! write(49, *) package(pack_index)%freq_cmf
   count_des_resd = count_des_resd + 1
   ldist = R_inf
   inCell = .FALSE.
   active = .FALSE.
   EXIT
   ! IF(chint < 1D-3) THEN
   !  ldist = norm2(package(pack_index)%pos - halfpos)
   !  ! write(*,*) 'resonance_distance: #1 ldist = ', ldist/R_inf
   !  ! STOP 'resonance_distance: testing'
   !  RETURN
   ! ELSE
   !  ! write(*,*) 'resonance_distance: pack_index = ', pack_index
   !  CALL abort()
   ! END IF
  END IF
  halfpos = upbond * 5.D-1 + lowbond * 5.D-1
  package(dummypackage)%pos = halfpos
  CALL doppler_factor(dummypackage, D)
  halffreq = package(pack_index)%freq_rf * D
  ! write(*,*) 'resonance_distance: lfreq / f = ', lfreq / f_line, &
  ! ' halffreq / f = ', halffreq / f_line, ' ufreq / f = ', ufreq / f_line
  ! write(*,*) 'resonance_distance: halffreq = ', halffreq
  ! decision which interval should we test next
  IF(upbond(1) == lowbond(1) .AND. upbond(2) == lowbond(2) &
   .AND. upbond(3) == lowbond(3)) THEN
   STOP 'resonance_distance: upbond == lowbond'
  END IF
  IF(halffreq > f_line) THEN
   chint = abs(halffreq - f_line) / f_line
   ! write(*,*) 'resonance_distance: chint = ', chint
   lowbond = halfpos
   lfreq = halffreq
   ! write(*,*) 'resonance_distance: lowbond = ', lowbond
  ELSE
   chint = abs(halffreq - f_line) / f_line
   ! write(*,*) 'resonance_distance: chint = ', chint
   upbond = halfpos
   ufreq = halffreq
   ! write(*,*) 'resonance_distance: upbond = ', upbond
  END IF
  IF(halffreq < bfreq .OR. halffreq > lfreq) THEN
   STOP 'resonance_distance: halffreq < bfreq or halffreq > lfreq'
  END IF
  ! write(*,*) 'resonance_distance: chint = ', chint, ' hf/bf = ', halffreq/bfreq,&
  ! '||left - right|| = ', norm2(lowbond - upbond), &
  ! '||pos - half|| = ', norm2(package(pack_index)%pos - halfpos)
  IF( chint <= minint ) THEN
   active = .FALSE.
   ldist = norm2(package(pack_index)%pos - halfpos)
   ! write(*,*) 'resonance_distance: #2 ldist = ', ldist/R_inf
  END IF
  IF(lfreq < f_line .OR. ufreq > f_line) STOP 'resonance_distance: freq do not fit'
 END DO ! until the calculation is active
 ! write(99,*) 'resonance_distance: I = ', I
 ! write(99,*) 'resonance_distance: ldista = ', light_speed * ( R_inf / V_inf ) * &
 !  ( ( package(pack_index)%freq_cmf / package(pack_index)%freq_rf)-&
 !   f_line / package(pack_index)%freq_rf) / R_star, ' cell_dist = ', cell_dist/R_star
 ! write(99,*) 'resonance_distance: ldist = ', ldist/R_star
! if the resonance point is not located in the actual propagation cell
! we have to go along the package path to the next cross boundary to
! realize if the resonant point is in this neighboor cell
ELSE
 inCell = .FALSE.
 ldist = R_inf
END IF

! write(*,*) 'resonance_distance: ldista = ', light_speed * ( R_inf / V_inf ) * &
!  ( ( package(pack_index)%freq_cmf / package(pack_index)%freq_rf)-&
!   f_line / package(pack_index)%freq_rf) / R_star, ' cell_dist = ', cell_dist/R_star
! write(*,*) 'resonance_distance: ldist = ', ldist/R_star

END SUBROUTINE resonance_distance
