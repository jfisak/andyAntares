SUBROUTINE next_cell_up(pack_index, dist, cell_down, next_cell)

USE types

IMPLICIT NONE

! input variables
DOUBLE PRECISION                        :: dist
INTEGER                                 :: cell_down, pack_index
! outpu variables
INTEGER                                 :: next_cell
! local variable next cell
! cross the surface
DOUBLE PRECISION, DIMENSION(3)          :: cross_pos
INTEGER                                 :: cross
! actual cell
INTEGER                                 :: act_cell
DOUBLE PRECISION, DIMENSION(3)          :: corner, width
INTEGER                                 :: upper_cell

act_cell = cell_down
cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist
cross = package(pack_index)%next_cross
IF(act_cell < 0) THEN
 next_cell = act_cell
 RETURN
END IF

print*, 'next_cell_up: act_cell = ', act_cell, ' cross_pos = ', cross_pos, 'cross = ', cross
DO
 corner = dyn_cell(act_cell)%corner
 width = dyn_cell(act_cell)%width
 upper_cell = dyn_cell(act_cell)%up_cell
 IF(dyn_cell(act_cell)%up_cell == 0) THEN
  next_cell = act_cell
  EXIT
 ELSE IF(dyn_cell(act_cell)%up_cell > 0) THEN
  ! we have to find which cell in the higher level corresponds to the cross point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in x direction
  IF(cross == posx .OR. cross == negx) THEN
   ! lower front cell
   IF(cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0 .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
    if(cross == posx) act_cell = upper_cell
    if(cross == negx) act_cell = upper_cell + 1
   ! lower rear cell
   ELSE IF(cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2) .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
    if(cross == posx) act_cell = upper_cell + 2
    if(cross == negx) act_cell = upper_cell + 3
   ! upper front cell
   ELSE IF(cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0 .AND. &
      cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
    if(cross == posx) act_cell = upper_cell + 4
    if(cross == posx) act_cell = upper_cell + 5
   ! upper rear cell
   ELSE IF(cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2) .AND. &
      cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
    if(cross == posx) act_cell = upper_cell + 6
    if(cross == posx) act_cell = upper_cell + 7
   END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in y direction
  ELSE IF(cross == posy .OR. cross == negy) THEN
   ! lower left cell
   IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
    if(cross == posy) act_cell = upper_cell
    if(cross == negy) act_cell = upper_cell + 2
   ! lower right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
    if(cross == posy) act_cell = upper_cell + 1
    if(cross == negy) act_cell = upper_cell + 3
   ! upper left cell
   ELSE IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
    if(cross == posy) act_cell = upper_cell + 4
    if(cross == posy) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
    if(cross == posy) act_cell = upper_cell + 5
    if(cross == posy) act_cell = upper_cell + 7
   END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in z direction
  ELSE IF(cross == posz .OR. cross == negz) THEN
   ! front left cell
   IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0) THEN
    if(cross == posz) act_cell = upper_cell
    if(cross == negz) act_cell = upper_cell + 4
   ! lower right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0) THEN
    if(cross == posz) act_cell = upper_cell + 1
    if(cross == negz) act_cell = upper_cell + 5
   ! upper left cell
   ELSE IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2)) THEN
    if(cross == posz) act_cell = upper_cell + 2
    if(cross == posz) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2)) THEN
    if(cross == posz) act_cell = upper_cell + 3
    if(cross == posz) act_cell = upper_cell + 7
   END IF
  END IF
 END IF  
END DO

END SUBROUTINE next_cell_up
