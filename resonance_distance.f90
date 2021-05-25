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
TYPE(photon)                    :: dummypackage
INTEGER                         :: dummy_pack
DOUBLE PRECISION                :: bfreq, halffreq, D
DOUBLE PRECISION                :: ufreq, lfreq
LOGICAL                         :: TESTING = .false.
LOGICAL                         :: iteration
LOGICAL                         :: inCell
LOGICAL                         :: isLdist
LOGICAL                         :: calculate
INTEGER                         :: I
INTEGER                         :: cell_number
DOUBLE PRECISION                :: dist, dist1, dist2
DOUBLE PRECISION                :: nu1, nu2
LOGICAL                         :: endit = .false.

INTEGER                         :: n_cell, next_cell, next_mi

DOUBLE PRECISION                :: ainx, binx, junk

inCell = .FALSE.
! IF(package(pack_index)%freq_cmf <= linelist(next1line)%freq) THEN
!  write(*,*) 'pack_index = ', pack_index
!  write(*,*) ' f_cmf / f_line = ', package(pack_index)%freq_cmf / linelist(next1line)%freq
!  STOP 'next_line'
! END IF

! write(78,*) pack_index, package(pack_index)%freq_rf, package(pack_index)%freq_cmf, &
!  package(pack_index)%freq_cmf / linelist(nextLine)%freq
! write(*,*) 'resonance_distance: pack_index = ', pack_index, ' f_line = ', f_line,&
!  ' cell_dist = ', cell_dist

minint = 1.D0 / linelist(1)%freq
calculate = .TRUE.
iteration = .TRUE.
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
 dummypackage = package(pack_index)
 package(pack_index)%pos(:) = package(pack_index)%pos(:) + ldist_analyt * package(pack_index)%dir(:)
 IF(ldist_analyt <= cell_dist) THEN
  inCell = .TRUE.
 ELSE
  inCell = .FALSE.
 END IF
 ldist = ldist_analyt
 RETURN
 ! write(*,*) 'resonance_distance: ldist_analyt = ', ldist_analyt / R_inf
END IF
! boundary coordinate
dummypackage = package(pack_index)
dummy_pack = SIZE(package)
package(dummy_pack) = dummypackage
CALL boundary3(dummy_pack, dist1, dist2, next_cell)
CALL move_package(dummy_pack, dist1)
package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D
package(pack_index)%e_cmf = package(pack_index)%e_rf * D
bfreq = dummypackage%freq_cmf
IF(velApprox /= 4) THEN
 rbound = dummypackage%pos
 ! frequency in the propagation cell boundary
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
    IF(I == 1000) THEN
     ! package(pack_index)%active = 0
     ! write(49, *) package(pack_index)%freq_cmf
     count_des_resd = count_des_resd + 1
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
    dummypackage%pos = halfpos
    CALL doppler_factor(0, dummypackage%pos, dummypackage%dir, D)
    halffreq = package(pack_index)%freq_rf * D
    ! IF(pack_index == 549) write(*,*) 'resonance_distance: lfreq / f = ', lfreq / f_line, &
    !  ' halffreq / f = ', halffreq / f_line, ' ufreq / f = ', ufreq / f_line
    ! write(*,*) 'resonance_distance: upbond = ', norm2(upbond), ' lowbond = ', norm2(lowbond)
    ! write(*,*) 'resonance_distance: halffreq = ', halffreq
    ! decision which interval should we test next
    IF(upbond(1) == lowbond(1) .AND. upbond(2) == lowbond(2) &
     .AND. upbond(3) == lowbond(3)) THEN
     ! write(78,*) 1.D8 * light_speed / package(pack_index)%freq_cmf
     write(*,*) 'resonance_distance: package = ', pack_index
     write(*,*) package(pack_index)%freq_cmf / linelist(nextLine)%freq
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
     CALL doppler_factor(0, dummypackage%pos, dummypackage%dir, D)
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
   dummypackage%pos = rbound
   lowbond = rbound
   CALL boundary3(dummypackage, dist, junk, cell_number)
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
   rbound = dummypackage%pos
   upbond = rbound
  END IF
 END DO
ELSE IF (velApprox == 4) THEN
 ! calculation of two boundary frequencies
 IF(next_cell > 0) THEN
  dummypackage%cell_numb = next_cell
  package(dummy_pack) = dummypackage
  CALL doppler_factor(dummy_pack, dummypackage%pos, dummypackage%dir, D)
  dummypackage%freq_cmf = dummypackage%freq_rf * D
  
  nu1 = package(pack_index)%freq_cmf
  nu2 = dummypackage%freq_cmf
   
  f_line = linelist(nextLine)%freq

  IF(nu1 < nu2) THEN
   inCell = .FALSE.
   RETURN
  END IF

  IF(nu1 > nu2) THEN
   IF(f_line > nu2 .AND. f_line < nu1) inCell = .TRUE.
  ELSE IF(nu1 < nu2) THEN
   IF(f_line < nu2 .AND. f_line > nu1) inCell = .TRUE.
  ELSE
   inCell = .FALSE.
  END IF

  if (inCell .EQV. .TRUE.) THEN
   ainx = (dist1 + dist2)/(nu2 - nu1)
   binx = (dist1 * nu1 + dist2 * nu2) / (nu1 - nu2)
   ldist = ainx * f_line + binx
   ! ldist = dist1
   ! write(*,*) 'resonance_distance: ainx = ', ainx, ' binx = ', binx/f_line
   ! write(*,*) 'resonance_distance: ldist = ', ldist/R_star, ' dist1 = ', dist1/R_star, ' inCell = ', inCell
   if (ldist < 0) then
    inCell = .FALSE.
   end if
  end if
 ELSE
  inCell = .FALSE.
 END IF

 
END IF

! write(*,*) 'resonance_distance: la = ', light_speed * ( R_inf / V_inf ) * &
!  ( ( package(pack_index)%freq_cmf / package(pack_index)%freq_rf)-&
!  f_line / package(pack_index)%freq_rf)/R_star, ' ldist = ', ldist/R_star

END SUBROUTINE resonance_distance
