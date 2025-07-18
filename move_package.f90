! Move photon package from the curent position for some distance (update package(pack_index)%pos)
!
! INPUT: pack_index(INT): index of a packet
!        dist(DBLE): distance to move
!        next_cell(INT): index of the next propGrid cell
!        change(LOG): should the index of the propGrid cell be changed?
! OUTPUT: NONE
!
SUBROUTINE move_package(pack_index, dist, next_cell, change)


USE types
USE constants
USE counters
USE dummypacket

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: dist, doppler_D

INTEGER                           :: next_cell
INTEGER                           :: cur_dummy_index
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: new_pos

LOGICAL                           :: change
DOUBLE PRECISION, PARAMETER       :: smallNumber = 1.D-2, mininum = 1.D4
! DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_corner, new_corner
! INTEGER                                                 :: cur_pgi
INTEGER                                                 :: next_cross
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: corner, upcorner, cur_pos
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: diff1, diff2
INTEGER, PARAMETER                                      :: corner_x = 1, corner_y = 2, corner_z = 3
INTEGER, PARAMETER                                      :: upcorner_x = 4, upcorner_y = 5, upcorner_z = 6
INTEGER                                                 :: cur_close
DOUBLE PRECISION                                        :: maxdist

IF(debug == 2) THEN
 write(*,*) 'move_package: going to move the packet, pack_index = ', pack_index
 write(*,*) 'move_package: dist = ', dist, ' dist/R_star = ', dist/R_star, 'dist/w = ', dist/basic_cell_width(ind_x)
END IF

IF(pack_index > SIZE(package)) THEN
 cur_dummy_index = pack_index - SIZE(package)
 dummypackage(cur_dummy_index)%pos(:) = dummypackage(cur_dummy_index)%pos(:) + dist * dummypackage(cur_dummy_index)%dir(:) 
 new_pos = dummypackage(cur_dummy_index)%pos
 if(abs(dummypackage(cur_dummy_index)%pos(ind_x)) < smallNumber) dummypackage(cur_dummy_index)%pos(ind_x) = 0.D0
 if(abs(dummypackage(cur_dummy_index)%pos(ind_y)) < smallNumber) dummypackage(cur_dummy_index)%pos(ind_y) = 0.D0
 if(abs(dummypackage(cur_dummy_index)%pos(ind_z)) < smallNumber) dummypackage(cur_dummy_index)%pos(ind_z) = 0.D0
ELSE
 package(pack_index)%pos(:) = package(pack_index)%pos(:) + dist * package(pack_index)%dir(:) 
 new_pos = package(pack_index)%pos
 next_cross = package(pack_index)%next_cross
 if(abs(package(pack_index)%pos(ind_x)) < smallNumber) package(pack_index)%pos(ind_x) = 0.D0
 if(abs(package(pack_index)%pos(ind_y)) < smallNumber) package(pack_index)%pos(ind_y) = 0.D0
 if(abs(package(pack_index)%pos(ind_z)) < smallNumber) package(pack_index)%pos(ind_z) = 0.D0
END IF

IF(debug == 2) THEN
 write(*,*) 'move_package: the packet was moved'
END IF



! Deactivate packets which travel beyond the photosphere
IF ((NORM2(new_pos) < R_star) .AND. (pack_index < SIZE(package))) THEN
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
! the change of packet is in this sbr available because during calling of the
! sbr doppler_factor a velocity vector is calculated and it should be done for
! the new propGrid cell index (if it was done for the old propGrid cell,
! it will generate a mistake)
IF(change) THEN
! we expect the packet to be exactly on the boundary, hence if the packet is only close to it
! we will slightly move it to be exactly there
 cur_pos = package(pack_index)%pos
 corner = dyn_cell(next_cell)%corner
 upcorner = dyn_cell(next_cell)%upcorner
 diff1 = abs(cur_pos - corner)
 diff2 = abs(cur_pos - upcorner)

 maxdist = 1.D99
 IF(diff1(ind_x) < maxdist) THEN
  maxdist = diff1(ind_x)
  cur_close = corner_x
 END IF
 IF(diff1(ind_y) < maxdist) THEN
  maxdist = diff1(ind_y)
  cur_close = corner_y
 END IF
 IF(diff1(ind_z) < maxdist) THEN
  maxdist = diff1(ind_z)
  cur_close = corner_z
 END IF
 IF(diff2(ind_x) < maxdist) THEN
  maxdist = diff2(ind_x)
  cur_close = upcorner_x
 END IF
 IF(diff2(ind_y) < maxdist) THEN
  maxdist = diff2(ind_y)
  cur_close = upcorner_y
 END IF
 IF(diff2(ind_z) < maxdist) THEN
  maxdist = diff2(ind_z)
  cur_close = upcorner_z
 END IF

 ! setting up a new coordinates
 IF(cur_close == corner_x) THEN
  package(pack_index)%pos(ind_x) = corner(ind_x)
 ELSE IF(cur_close == corner_y) THEN
  package(pack_index)%pos(ind_y) = corner(ind_y)
 ELSE IF(cur_close == corner_z) THEN
  package(pack_index)%pos(ind_z) = corner(ind_z)
 ELSE IF(cur_close == upcorner_x) THEN
  package(pack_index)%pos(ind_x) = upcorner(ind_x)
 ELSE IF(cur_close == upcorner_y) THEN
  package(pack_index)%pos(ind_y) = upcorner(ind_y)
 ELSE IF(cur_close == upcorner_z) THEN
  package(pack_index)%pos(ind_z) = upcorner(ind_z)
 END IF

 
 CALL change_cell(pack_index, next_cell)
END IF
CALL doppler_factor(pack_index, doppler_D)
IF(pack_index <= SIZE(package)) THEN
 package(pack_index)%freq_cmf = package(pack_index)%freq_rf * doppler_D
 package(pack_index)%e_cmf = package(pack_index)%e_rf * doppler_D
 package(pack_index)%delta_s = package(pack_index)%delta_s + dist
ELSE IF(pack_index > SIZE(package)) THEN
 dummypackage(cur_dummy_index)%freq_cmf = dummypackage(cur_dummy_index)%freq_rf * doppler_D
 dummypackage(cur_dummy_index)%e_cmf = dummypackage(cur_dummy_index)%e_rf * doppler_D
END IF
IF(pack_index <= SIZE(package)) THEN
 IF (package(pack_index)%freq_cmf < 0) THEN
  write(*,*) 'move_package: package = ', pack_index, ' prop. cell = ', package(pack_index)%cell_numb
  write(*,*) 'FREQUENCY IS LOWER THAN ZERO!!!'
 END IF
END IF


END SUBROUTINE move_package
