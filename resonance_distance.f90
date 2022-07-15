SUBROUTINE resonance_distance(pack_index, nextLine, f_line, cell_dist, ldist, inCell, isLdist)
USE types
USE counters
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index, nextLine
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
INTEGER                         :: dummypackage
DOUBLE PRECISION                :: bfreq, halffreq, D
DOUBLE PRECISION                :: ufreq, lfreq
LOGICAL                         :: TESTING = .false.
LOGICAL                         :: iteration
LOGICAL                         :: inCell
LOGICAL                         :: isLdist
LOGICAL                         :: calculate
INTEGER                         :: I
INTEGER                         :: cell_number
DOUBLE PRECISION                :: dist
LOGICAL                         :: endit = .false.

! IF(package(pack_index)%freq_cmf <= linelist(next1line)%freq) THEN
!  write(*,*) 'pack_index = ', pack_index
!  write(*,*) ' f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq
!  STOP 'next_line'
! END IF

! write(78,*) pack_index, package(pack_index)%freq_rf, package(pack_index)%freq_cmf, &
!  package(pack_index)%freq_cmf / linelist(nextLine)%freq
dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
! write(*,*) 'resonance_distance: pack_index = ', pack_index, ' f_line = ', f_line,&
!  ' cell_dist = ', cell_dist

! write(*,*) '*****************************************************'
! 
! IF(pack_index == 549) write(*,*) 'resonance_distance: fcmf / f_line = ', package(pack_index)%freq_cmf / f_line

! write(*,*) 'resonance_distance***************************************'
! minint = 1.D0 / linelist(1)%freq
calculate = .TRUE.
iteration = .TRUE.
minint = 1.D-4
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
 CALL move_package(dummypackage, ldist_analyt)
 IF(ldist_analyt <= cell_dist) THEN
  inCell = .TRUE.
 ELSE
  inCell = .FALSE.
 END IF
 ldist = ldist_analyt
 RETURN
 ! write(*,*) 'resonance_distance: ldist_analyt = ', ldist_analyt / R_inf
END IF

IF(velapprox == 2) THEN
 CALL move_package(dummypackage, cell_dist)
 bfreq = package(dummypackage)%freq_cmf
 IF(package(pack_index)%freq_cmf > f_line .AND. bfreq < f_line) THEN
  inCell = .TRUE.
  ldist = cell_dist/2.0
 ELSE
  inCell = .FALSE.
 END IF
 RETURN
END IF


 ! boundary coordinate
! write(*,*) 'resonance_distance: cell_dist = ', cell_dist
CALL move_package(dummypackage, cell_dist)
bfreq = package(dummypackage)%freq_cmf
rbound = package(dummypackage)%pos
! frequency in the propagation cell boundary
! package(dummypackage)%pos = rbound
lowbond = package(pack_index)%pos
DO WHILE(iteration)
 lfreq = package(pack_index)%freq_cmf
 upbond = rbound
 ! cmf frequency on boundary
 IF(lfreq < f_line) THEN
  calculate = .FALSE.
  iteration = .FALSE.
 END IF
 ufreq = bfreq
 ! we can calculate the line_dist
 ! lower and upper boundary
 ! IF(pack_index == 549) write(*,*) 'resonance_distance: bfreq / f_line = ', bfreq / f_line
 IF(bfreq < f_line .AND. calculate .EQV. .TRUE.) THEN
  active = .TRUE.
  inCell = .TRUE.
  ! write(*,*) 'resonance_distance: line res in prop cell'
  I = 0
  DO WHILE(active .EQV. .TRUE.)
   I = I + 1
   ! write(*,*) 'resonance_distance: I = ', I
   IF(I == 100) THEN
    ! package(pack_index)%active = 0
    ! write(49, *) package(pack_index)%freq_cmf
    ! write(*,*) 'resonance_distance: too many iterations...'
    ! write(*,*) 'resonance_distance: packet: ', pack_index, ' too many iterations...'
    ! STOP    
    ! count_des_resd = count_des_resd + 1
    ldist = R_inf
    inCell = .FALSE.
    active = .FALSE.
    iteration = .FALSE.
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
   ! IF(pack_index == 549)
   ! write(*,*) 'resonance_distance: lfreq / f = ', lfreq / f_line
   ! write(*,*)  ' halffreq / f = ', halffreq / f_line, ' ufreq / f = ', ufreq / f_line
   ! write(*,*) 'resonance_distance: upbond = ', norm2(upbond), ' lowbond = ', norm2(lowbond)
   ! write(*,*) 'resonance_distance: halffreq = ', halffreq
   ! decision which interval should we test next
   IF(upbond(1) == lowbond(1) .AND. upbond(2) == lowbond(2) &
    .AND. upbond(3) == lowbond(3)) THEN
    ! write(78,*) 1.D8 * light_speed / package(pack_index)%freq_cmf
    ! write(*,*) 'resonance_distance: package = ', pack_index
    ! write(*,*) package(pack_index)%freq_cmf / linelist(nextLine)%freq
    STOP 'upbond == lowbond'
    active = .FALSE.
    iteration = .FALSE.
    ! inCell = .FALSE.
    ldist = norm2(package(pack_index)%pos - halfpos)
    ! ldist = R_inf
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
   IF(halffreq < bfreq) THEN
    write(*,*) 'resonance_distance: halffreq < bfreq'
    endit = .TRUE.
   END IF
   IF(halffreq > lfreq) THEN
    write(*,*) 'resonance_distance: halffreq > lfreq'
    endit = .TRUE.
   END IF
   IF(endit) THEN
    write(*,*) 'resonance_distance: pack_index = ', pack_index, ' f_line = ', f_line,&
     ' cell_dist = ', cell_dist/ R_inf
    write(*,*) 'resonance_distance: r / R_inf = ', norm2(lowbond) / R_inf
    write(*,*) 'resonance_distance: r_r / R_inf = ', norm2(upbond) / R_inf
    ! write(*,*) 'resonance_distance: r / R_star = ', norm2(halfpos) / R_star
    CALL doppler_factor(dummypackage, D)
    write(*,*) 'resonance_distance: D = ', D
    STOP 'resonance_distance'
   END IF
   IF( chint <= minint ) THEN
    active = .FALSE.
    ldist = norm2(package(pack_index)%pos - halfpos)
    ! write(*,*) 'resonance_distance: #2 ldist = ', ldist/R_inf
   END IF
   IF(lfreq < f_line .OR. ufreq > f_line) STOP 'resonance_distance: freq do not fit'
  END DO ! until the calculation is active
  iteration = .FALSE.
 ! if the resonance point is not located in the actual propagation cell
 ! we have to go along the package path to the next cross boundary to
 ! realize if the resonant point is in this neighboor cell
 ELSE
  !write(*,*) 'resonance_distance: line res out off prop cell'
  ! only continnum process could happen
  IF(isLdist) THEN
   ldist = R_inf
   inCell = .FALSE.
   RETURN
  END IF
  package(dummypackage)%pos = rbound
  lowbond = rbound
  CALL boundary3(dummypackage, dist, cell_number)
  IF(cell_number <= 0) THEN
   inCell = .FALSE.
   ldist = R_inf
   RETURN
  END IF
  CALL change_cell(dummypackage, cell_number)
  CALL move_package(dummypackage, dist)
  ! write(*,*) 'resonance_distance: rbound / ldist_analyt = ',&
  ! norm2(rbound) / ldist_analyt
  ! IF(norm2(rbound - package(pack_index)%pos) > ldist_analyt) write(*,*) 'resonance_distance: ||rbound|| > ldist_analyt'
  IF(norm2(rbound) > R_inf) THEN
   ldist = R_inf
   inCell = .FALSE.
   ! write(*,*) 'resonance_distance: #3 ldist = ', ldist/R_inf
   RETURN
  END IF
  rbound = package(dummypackage)%pos
  upbond = rbound
   ! write(*,*) 'resonance_distance: dist = ', norm2(package(pack_index)%pos - rbound)/R_inf, norm2(package(dummypackage)%pos)/R_inf
   ! write(*,*) 'resonance_distance: rbound = ', norm2(rbound)/R_inf
   ! STOP 'resonance_distance: testing'
 END IF
END DO

! write(*,*) 'resonance_distance: la = ', light_speed * ( R_inf / V_inf ) * &
!  ( ( package(pack_index)%freq_cmf / package(pack_index)%freq_rf)-&
!  f_line / package(pack_index)%freq_rf)/R_star, ' ldist = ', ldist/R_star

END SUBROUTINE resonance_distance
