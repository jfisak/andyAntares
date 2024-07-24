SUBROUTINE move_package(pack_index, dist, next_cell, change)

! Move photon package from the curent position for some distance (update package(pack_index)%pos)

USE types
USE constants
USE counters
USE dummypacket

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: dist, doppler_D
INTEGER                           :: ind_dummypackage

INTEGER                           :: next_cell
INTEGER                           :: cur_dummy_index
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: new_pos

LOGICAL                           :: change
DOUBLE PRECISION, PARAMETER       :: smallNumber = 1.D-2

ind_dummypackage = SIZE(package)

IF(pack_index > SIZE(package)) THEN
 cur_dummy_index = pack_index - SIZE(package)
 dummypackage(pack_index)%pos(:) = dummypackage(pack_index)%pos(:) + dist * dummypackage(pack_index)%dir(:) 
 new_pos = dummypackage(pack_index)%pos
 if(abs(dummypackage(pack_index)%pos(ind_x)) < smallNumber) dummypackage(pack_index)%pos(ind_x) = 0e0
 if(abs(dummypackage(pack_index)%pos(ind_y)) < smallNumber) dummypackage(pack_index)%pos(ind_y) = 0e0
 if(abs(dummypackage(pack_index)%pos(ind_z)) < smallNumber) dummypackage(pack_index)%pos(ind_z) = 0e0
ELSE
 package(pack_index)%pos(:) = package(pack_index)%pos(:) + dist * package(pack_index)%dir(:) 
 new_pos = package(pack_index)%pos
 if(abs(package(pack_index)%pos(ind_x)) < smallNumber) package(pack_index)%pos(ind_x) = 0e0
 if(abs(package(pack_index)%pos(ind_y)) < smallNumber) package(pack_index)%pos(ind_y) = 0e0
 if(abs(package(pack_index)%pos(ind_z)) < smallNumber) package(pack_index)%pos(ind_z) = 0e0
END IF


! Deactivate packets which travel beyond the photosphere
IF ((NORM2(new_pos) < R_star) .AND. (pack_index .NE. ind_dummypackage)) THEN
 IF(abs_surface >= 1) THEN
  CALL photosphere_interaction(pack_index)
  dist = 0.0
  change=.false.
 ELSE
  IF(pack_index <= SIZE(package)) THEN
   package(pack_index)%active = 0
  END IF
  count_des_phot = count_des_phot + 1
 END IF
END IF
! Rest frame quantities do not change while propagating without any events, 
! but cmf quantities need to be updated
IF(change) CALL change_cell(pack_index, next_cell)
CALL doppler_factor(pack_index, doppler_D)
IF(pack_index <= SIZE(package)) THEN
 package(pack_index)%freq_cmf = package(pack_index)%freq_rf * doppler_D
 package(pack_index)%e_cmf = package(pack_index)%e_rf * doppler_D
 package(pack_index)%delta_s = package(pack_index)%delta_s + dist
ELSE IF(pack_index > SIZE(package)) THEN
 dummypackage(pack_index)%freq_cmf = dummypackage(pack_index)%freq_rf * doppler_D
 dummypackage(pack_index)%e_cmf = dummypackage(pack_index)%e_rf * doppler_D
END IF
IF (package(pack_index)%freq_cmf < 0) THEN
 write(*,*) 'move_package: package = ', pack_index, ' prop. cell = ', package(pack_index)%cell_numb
 write(*,*) 'FREQUENCY IS LOWER THAN ZERO!!!'
END IF

END SUBROUTINE move_package
