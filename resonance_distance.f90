SUBROUTINE resonance_distance(pack_index, nextLine, cell_dist, ldist)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index, nextLine
DOUBLE PRECISION                :: cell_dist
DOUBLE PRECISION                :: f_line
! output variables
DOUBLE PRECISION                :: ldist
! beta law 
LOGICAL                         :: active = .TRUE.
DOUBLE PRECISION, DIMENSION(3)  :: halfpos, lowbond, upbond, rbound
! DOUBLE PRECISION, PARAMETER     :: minint = 1.D-15
DOUBLE PRECISION                :: minint
DOUBLE PRECISION                :: chint
! testing
DOUBLE PRECISION                :: ldist_analyt
INTEGER                         :: my_rank, n_pack_d, dummypackage
DOUBLE PRECISION                :: bfreq, halffreq, D
INTEGER                         :: OMP_GET_THREAD_NUM
DOUBLE PRECISION                :: ufreq, lfreq
LOGICAL                         :: TESTING = .TRUE.
LOGICAL                         :: outOfCell

my_rank = OMP_GET_THREAD_NUM()
n_pack_d = SIZE(package)
dummypackage = n_pack_d - n_dummy_packs + my_rank + 1
package(dummypackage) = package(pack_index)

f_line = linelist(nextLine)%freq

! write(*,*) 'resonance_distance: nextLine = ', nextLine

minint = 1.D-1 / linelist(1)%freq
! Calculate distance the photon needs to travel to come to
! resonance with the next line. This assumes homologous
! expansion i.e. velocity is proportional to r. Projected
! gradient of the projected velocity in a direction of the
! photon propagation is more complicated in case of no
! homologous expansion
! homologous approximation
! IF(velapprox == 0) THEN
!  ldist = light_speed * (R_inf/V_inf) * &
!   ((package(pack_index)%freq_cmf - f_line) / package(pack_index)%freq_rf)
! ELSE
 IF(TESTING .EQV. .TRUE.) THEN
  ldist_analyt = light_speed * ( R_inf / V_inf ) * &
   ( ( package(pack_index)%freq_cmf - f_line ) / package(pack_index)%freq_rf)
 END IF
  ! boundary coordinate
 ! write(*,*) 'resonance_distance: cell_dist = ', cell_dist
 rbound = package(pack_index)%pos + package(pack_index)%dir * cell_dist
 ! frequency in the propagation cell boundary
 package(dummypackage)%pos = rbound
 CALL doppler_factor(dummypackage, D)
 bfreq = package(pack_index)%freq_rf * D
 ! write(*,*) 'resonance_distance: rbound = ', rbound
 lfreq = package(pack_index)%freq_cmf
 ufreq = bfreq
 ! we can calculate the line_dist
 ! lower and upper boundary
 lowbond = package(pack_index)%pos
 upbond = rbound
 ! write(*,*) 'bfreq / f_line = ', bfreq / f_line
 IF(bfreq < f_line) THEN
  outOfCell = .FALSE.
  ! write(*,*) 'resonance_distance: line res in prop cell'
  DO WHILE(active .EQV. .TRUE.)
   ! what is the value in the half of the interval
   halfpos = (upbond + lowbond) / 2.D0
   ! write(*,*) 'resonance_distance: upbond = ', upbond, 'lowbond = ', lowbond, ' halfpos = ', halfpos
   package(dummypackage)%pos = halfpos
   CALL doppler_factor(dummypackage, D)
   halffreq = package(pack_index)%freq_rf * D
   ! decision which interval should we test next
   IF(halffreq > f_line) THEN
    chint = abs(halffreq - f_line) / f_line
    ! write(*,*) 'resonance_distance: chint = ', chint
    lowbond = halfpos
    ! write(*,*) 'resonance_distance: lowbond = ', lowbond
   ELSE
    chint = abs(halffreq - f_line) / f_line
    ! write(*,*) 'resonance_distance: chint = ', chint
    upbond = halfpos
    ! write(*,*) 'resonance_distance: upbond = ', upbond
   END IF
   IF( chint < minint ) THEN
    active = .FALSE.
    ldist = norm2(package(pack_index)%pos - halfpos)
   END IF
  END DO ! until the calculation is active
  IF(TESTING .EQV. .TRUE.) write(*,*) 'resonance_distance: ldist_analyt = ',&
   ldist_analyt, ' ldist = ', ldist, abs(1.0 - ldist_analyt / ldist)
  ! write(*,*) 'resonance_distance: f_line = ', f_line / 1.D14, &
  !  ' halffreq = ', halffreq/1.D14, ' fl / hf = ', abs(halffreq - f_line) / f_line
  ! STOP 'resonance_distance: testing'
 ELSE
  ! write(*,*) 'resonance_distance: line res out off prop cell'
  ! only continnum process could happen
  outOfCell = .TRUE.
  ldist = 1.D20
 END IF
 ! coordinates of a line interaction
! END IF

! calculation of the line optical depth


END SUBROUTINE resonance_distance
