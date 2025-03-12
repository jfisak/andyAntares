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

IF(debug == 2) THEN
 write(*,*) 'move_package: going to move the packet, pack_index = ', pack_index
 write(*,*) 'move_package: dist = ', dist
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

 ! correction of a position
 ! IF(next_cell > 0) THEN
 !  cur_pgi = package(pack_index)%cell_numb
 !  new_pos = package(pack_index)%pos
 !  cur_corner = dyn_cell(cur_pgi)%corner
 !  new_corner = dyn_cell(next_cell)%corner
 !  IF(next_cross == posx) THEN
 !   IF(new_pos(ind_x) > new_corner(ind_x) - mininum .and. new_pos(ind_x) < new_corner(ind_x) + mininum) THEN
 !    package(pack_index)%pos(ind_x) = new_corner(ind_x)
 !   END IF
 !  ELSE IF(next_cross == posy) THEN
 !   IF(new_pos(ind_y) > new_corner(ind_y) - mininum .and. new_pos(ind_y) < new_corner(ind_y) + mininum) THEN
 !    package(pack_index)%pos(ind_y) = new_corner(ind_y)
 !   END IF
 !  ELSE IF(next_cross == posz) THEN
 !   IF(new_pos(ind_z) > new_corner(ind_z) - mininum .and. new_pos(ind_z) < new_corner(ind_z) + mininum) THEN
 !    package(pack_index)%pos(ind_z) = new_corner(ind_z)
 !   END IF
 !  ELSE IF(next_cross == negx) THEN
 !   IF(new_pos(ind_x) > cur_corner(ind_x) - mininum .and. new_pos(ind_x) < cur_corner(ind_x) + mininum) THEN
 !    package(pack_index)%pos(ind_x) = cur_corner(ind_x)
 !   END IF
 !  ELSE IF(next_cross == negy) THEN
 !   IF(new_pos(ind_y) > cur_corner(ind_y) - mininum .and. new_pos(ind_y) < cur_corner(ind_y) + mininum) THEN
 !    package(pack_index)%pos(ind_y) = cur_corner(ind_y)
 !   END IF
 !  ELSE IF(next_cross == negz) THEN
 !   IF(new_pos(ind_z) > cur_corner(ind_z) - mininum .and. new_pos(ind_z) < cur_corner(ind_z) + mininum) THEN
 !    package(pack_index)%pos(ind_z) = cur_corner(ind_z)
 !   END IF
 !  END IF
 ! END IF ! next_cell > 0
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
IF(change) CALL change_cell(pack_index, next_cell)
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
