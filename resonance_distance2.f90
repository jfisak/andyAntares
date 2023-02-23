! this sbr calculates if a packet can interact with the next line in the current propagation cell
SUBROUTINE resonance_distance2(pack_index, nextLine, cell_dist, inCell, ldist)

USE types
USE constants

IMPLICIT NONE

INTEGER                                 :: pack_index, nextLine
DOUBLE PRECISION                        :: cell_dist
LOGICAL                                 :: inCell, redshift
DOUBLE PRECISION                        :: ldist


INTEGER                                 :: dummypackage, cell_number
LOGICAL                                 :: iteration, change_of_cell
DOUBLE PRECISION                        :: minint

DOUBLE PRECISION                        :: bfreq, lfreq, ufreq, f_line
DOUBLE PRECISION, DIMENSION(3)          :: rbond, lbond, ubond
DOUBLE PRECISION                        :: halffreq
DOUBLE PRECISION, DIMENSION(3)          :: halfpos

DOUBLE PRECISION                        :: D
DOUBLE PRECISION                        :: chint

INTEGER                                 :: I
INTEGER, PARAMETER                      :: maxit = 111
! test of convergence


redshift = package(pack_index)%redshift
f_line = linelist(nextLine)%freq

dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)
cell_number = package(pack_index)%cell_numb
minint = 1.D-3


! we just move dummypackage to the boundary to calculate its CMF frequency
! a change of cell is not required
change_of_cell = .FALSE.
CALL move_package(dummypackage, cell_dist, cell_number, change_of_cell)

! forward boundaries
bfreq = package(dummypackage)%freq_cmf
rbond = package(dummypackage)%pos
ufreq = bfreq
ubond = rbond

! current boundaries (in the current packet position)
lbond = package(pack_index)%pos
lfreq = package(pack_index)%freq_cmf

iteration = .TRUE.
inCell = .TRUE.


! test if the line frequency is in the interval
IF(redshift) THEN
 ! frequency should getting lower
 IF(f_line < bfreq .or. f_line > lfreq) THEN
  inCell = .false.
  ldist = R_inf
  RETURN
 END IF
ELSE ! blueshift
 IF(f_line > bfreq .or. f_line < lfreq) THEN
  inCell = .false.
  ldist = R_inf
  RETURN
 END IF
END IF



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
 package(dummypackage)%pos = halfpos
 CALL doppler_factor(dummypackage, D)
 halffreq = package(pack_index)%freq_rf * D

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



























END SUBROUTINE resonance_distance2

