! this sbr calculates if a packet can interact with the next line in the current propagation cell
SUBROUTINE resonance_distance2(pack_index, nextLine, cell_dist, inCell, ldist)

USE types
USE constants
USE dummypacket

IMPLICIT NONE

INTEGER                                 :: pack_index, nextLine
INTEGER                                 :: cur_dummypack, dummypack_index
DOUBLE PRECISION                        :: cell_dist
LOGICAL                                 :: inCell, redshift
DOUBLE PRECISION                        :: ldist


INTEGER                                 :: cell_number
LOGICAL                                 :: iteration, change_of_cell
DOUBLE PRECISION, PARAMETER             :: minint = 1.D-3

DOUBLE PRECISION                        :: bfreq, lfreq, ufreq, f_line
DOUBLE PRECISION, DIMENSION(3)          :: rbond, lbond, ubond
DOUBLE PRECISION                        :: halffreq
DOUBLE PRECISION, DIMENSION(3)          :: halfpos

DOUBLE PRECISION                        :: D_doppler
DOUBLE PRECISION                        :: chint

INTEGER                                 :: I
INTEGER, PARAMETER                      :: maxit = 111
! test of convergence


redshift = package(pack_index)%redshift
f_line = linelist(nextLine)%freq
! write(*,*) 'resonance_distance2: redshift = ', redshift

cur_dummypack = find_free_index() 
CALL copy_package(pack_index, cur_dummypack)
dummypack_index = cur_dummypack + SIZE(package)

cell_number = package(pack_index)%cell_numb

! we just move dummypackage to the boundary to calculate its CMF frequency
! a change of cell is not required
change_of_cell = .FALSE.
CALL move_package(dummypack_index, cell_dist, cell_number, change_of_cell)

! forward boundaries
bfreq = dummypackage(cur_dummypack)%freq_cmf
! write(*,*) 'resonance_distance2: f_rf = ', package(dummypackage)%freq_rf, 'f_cmf/f_line = ', package(dummypackage)%freq_cmf/f_line
rbond = dummypackage(cur_dummypack)%pos
ufreq = bfreq
ubond = rbond

! **testing**
! ldist = const_c * (R_inf/V_inf) * ((package(pack_index)%freq_cmf - f_line)/f_line)
! if(ldist > cell_dist) then
!  ldist = R_inf
! end if
! return
! **testing**

! current boundaries (in the current packet position)
lbond = package(pack_index)%pos
lfreq = package(pack_index)%freq_cmf

iteration = .TRUE.
inCell = .TRUE.

! write(*,*) 'resonance_distance2: f_line/f_cmf = ', f_line/ufreq
! test if the line frequency is in the interval
IF(redshift) THEN
 ! frequency should getting lower
 IF(f_line < ufreq .or. f_line > lfreq) THEN
  ! write(*,*) 'resonance_distance2: packet = ', pack_index, ' f_line < bfreq .or. f_line > lfreq'
  inCell = .false.
  ldist = R_inf
  CALL deactivate_dummy_packet(cur_dummypack)
  RETURN
 END IF
ELSE ! blueshift
 IF(f_line > bfreq .or. f_line < lfreq) THEN
  inCell = .false.
  ldist = R_inf
  CALL deactivate_dummy_packet(cur_dummypack)
  RETURN
 END IF
END IF


! write(*,*) 'resonance_distance2: calculation of a resonance point'
! write(*,*) 'resonance_distance2: inCell = ', inCell
I = 0
DO WHILE(iteration)
 
 I = I + 1
 if(I == maxit) then
  ldist = norm2(package(pack_index)%pos - halfpos)
  iteration = .false.
  ! STOP 'resonance_distance2: testing'
 end if


 ! a decision which interval will be hit -- packets or boundary
 ! moving a virtual packet to the center between packet and boundary
 halfpos = ubond * 5.D-1 + lbond * 5.D-1
 ! write(*,*) 'resonance_distance2: halfpos = ', halfpos
 dummypackage(cur_dummypack)%pos = halfpos
 CALL doppler_factor(dummypack_index, D_doppler)
 halffreq = package(pack_index)%freq_rf * D_doppler

 ! test of convergence
 chint = abs(halffreq - f_line) / f_line
 IF(chint <= minint) THEN
  iteration = .false.
  ldist = norm2(package(pack_index)%pos - halfpos)
 END IF


 ! the decision itself
 ! depends on the packet red or blue shift
 IF(redshift) THEN
  if(f_line > halffreq) then
   ubond = halfpos
   ufreq = halffreq
  else
   lbond = halfpos
   lfreq = halffreq
  end if
 ELSE ! blueshift
  if(f_line > halffreq) then
   lbond = halfpos
   lfreq = halffreq
  else
   ubond = halfpos
   ufreq = halffreq
  end if
 END IF





END DO

CALL deactivate_dummy_packet(cur_dummypack)
! write(*,*) 'resonance_distance2: after calc, inCell = ', inCell

END SUBROUTINE resonance_distance2

