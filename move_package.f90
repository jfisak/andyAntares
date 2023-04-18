SUBROUTINE move_package(pack_index, dist, next_cell, change)

! Move photon package from the curent position for some distance (update package(pack_index)%pos)

USE types
USE constants
USE counters

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: dist, D, vec_length
INTEGER                           :: dummypackage

INTEGER                           :: next_cell

LOGICAL                           :: change

dummypackage = SIZE(package)

! write(*,*) 'move_package: pack_index = ', pack_index, ' dist = ', dist/R_star

! Calculate the position of package
package(pack_index)%pos(1) = package(pack_index)%pos(1) + dist * package(pack_index)%dir(1) 
package(pack_index)%pos(2) = package(pack_index)%pos(2) + dist * package(pack_index)%dir(2) 
package(pack_index)%pos(3) = package(pack_index)%pos(3) + dist * package(pack_index)%dir(3)
if(abs(package(pack_index)%pos(1)) < 1e-1) package(pack_index)%pos(1) = 0e0
if(abs(package(pack_index)%pos(2)) < 1e-1) package(pack_index)%pos(2) = 0e0
if(abs(package(pack_index)%pos(3)) < 1e-1) package(pack_index)%pos(3) = 0e0

! Deactivate packets which travel beyond the photosphere
IF ((vec_length(package(pack_index)%pos) < R_star) .AND. (pack_index .NE. dummypackage)) THEN
 IF(abs_surface >= 1) THEN
  CALL photosphere_interaction(pack_index)
  dist = 0.0
  change=.false.
 ELSE
  package(pack_index)%active = 0
  count_des_phot = count_des_phot + 1
 END IF
END IF
! Rest frame quantities do not change while propagating without any events, 
! but cmf quantities need to be updated
IF(change) CALL change_cell(pack_index, next_cell)
CALL doppler_factor(pack_index, D)
package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D
package(pack_index)%e_cmf = package(pack_index)%e_rf * D
IF (package(pack_index)%freq_cmf < 0) THEN
 write(*,*) 'move_package: package = ', pack_index, ' prop. cell = ', package(pack_index)%cell_numb
 write(*,*) 'FREQUENCY IS LOWER THAN ZERO!!!'
 !STOP 
END IF
package(pack_index)%delta_s = package(pack_index)%delta_s + dist
! write(79,*) package(pack_index)%delta_s/R_star, package(pack_index)%freq_cmf

END SUBROUTINE move_package
