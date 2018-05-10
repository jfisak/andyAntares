SUBROUTINE resonance_distance(pack_index, nextLine, cell_dist, ldist, inCell)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index, nextLine
DOUBLE PRECISION                :: cell_dist
DOUBLE PRECISION                :: f_line
! output variables
DOUBLE PRECISION                :: ldist
! beta law 
LOGICAL                         :: active
DOUBLE PRECISION, DIMENSION(3)  :: halfpos, lowbond, upbond, rbound, lobound
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
LOGICAL                         :: outOfCell
LOGICAL                         :: inCell
INTEGER                         :: I
INTEGER                         :: cell_number
DOUBLE PRECISION                :: dist

dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
outOfCell = .TRUE.

f_line = linelist(nextLine)%freq
! write(*,*) '*****************************************************'

IF(package(pack_index)%freq_cmf < f_line) THEN
 ldist = 0.D0
 RETURN
END IF

! write(*,*) 'resonance_distance***************************************'
! minint = 1.D0 / linelist(1)%freq
minint = 1.D-5
! Calculate distance the photon needs to travel to come to
! resonance with the next line. This assumes homologous
! expansion i.e. velocity is proportional to r. Projected
! gradient of the projected velocity in a direction of the
! photon propagation is more complicated in case of no
! homologous expansion
! homologous approximation
 IF(TESTING .EQV. .TRUE.) THEN
  ldist_analyt = light_speed * ( R_inf / V_inf ) * &
   ( ( package(pack_index)%freq_cmf - f_line ) / package(pack_index)%freq_rf)
  ! write(*,*) 'resonance_distance: ldist_analyt = ', ldist_analyt / R_inf
 END IF
  ! boundary coordinate
 ! write(*,*) 'resonance_distance: cell_dist = ', cell_dist
 lobound = package(pack_index)%pos
 rbound = package(pack_index)%pos + package(pack_index)%dir * cell_dist
 ! frequency in the propagation cell boundary
 package(dummypackage)%pos = rbound
 lowbond = package(pack_index)%pos
 upbond = rbound
 ! write(*,*) 'resonance_distance: outOfCell = ', outOfCell
DO WHILE(outOfCell .EQV. .TRUE.)
 ! cmf frequency on boundary
 CALL doppler_factor(dummypackage, D)
 bfreq = package(dummypackage)%freq_rf * D
 lfreq = package(pack_index)%freq_cmf
 IF(lfreq < f_line) STOP 'resonance_distance: lfreq < f_line'
 ufreq = bfreq
 ! we can calculate the line_dist
 ! lower and upper boundary
 ! write(*,*) 'resonance_distance: bfreq / f_line = ', bfreq / f_line
 IF(bfreq < f_line) THEN
  active = .TRUE.
  outOfCell = .FALSE.
  inCell = .TRUE.
  ! write(*,*) 'resonance_distance: line res in prop cell'
  I = 0
  DO WHILE(active .EQV. .TRUE.)
   I = I + 1
   ! write(*,*) 'resonance_distance: I = ', I
   IF(I == 1000) THEN
    IF(chint < 1D1 * minint) THEN
     ldist = norm2(package(pack_index)%pos - halfpos)
     ! write(*,*) 'resonance_distance: #1 ldist = ', ldist/R_inf
     RETURN
    ELSE
     CALL abort()
    END IF
   END IF
   halfpos = upbond * 5.D-1 + lowbond * 5.D-1
   package(dummypackage)%pos = halfpos
   CALL doppler_factor(dummypackage, D)
   halffreq = package(pack_index)%freq_rf * D
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
   IF( chint <= minint ) THEN
    active = .FALSE.
    ldist = norm2(package(pack_index)%pos - halfpos)
    ! write(*,*) 'resonance_distance: #2 ldist = ', ldist/R_inf
   END IF
   IF(lfreq < f_line .OR. ufreq > f_line) STOP 'resonance_distance: freq do not fit'
  END DO ! until the calculation is active
 ! if the resonance point is not located in the actual propagation cell
 ! we have to go along the package path to the next cross boundary to
 ! realize if the resonant point is in this neighboor cell
 ELSE
  !write(*,*) 'resonance_distance: line res out off prop cell'
  ! only continnum process could happen
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
END DO ! while it is not located in cell

! write(*,*) 'resonance_distance: #4 ldist = ', ldist/R_inf
END SUBROUTINE resonance_distance
