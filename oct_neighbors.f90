SUBROUTINE oct_neighbors(pack_index, rel_pos, velgridcells)


USE types
IMPLICIT NONE

INTEGER, PARAMETER                      :: n_oct=8

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: rel_pos

INTEGER                                 :: act_cell
DOUBLE PRECISION, DIMENSION(3)          :: act_corner, act_width

DOUBLE PRECISION, DIMENSION(3)          :: point_width, point_corner, point_vec
INTEGER                                 :: point_cell, cur_cell

INTEGER, DIMENSION(n_oct)                   :: velgridcells

INTEGER                                 :: I, J, dummy, a


act_cell = package(pack_index)%cell_numb
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width

! this part code only for the regular grid
IF(dyngrid == 0) THEN
 ! find the neighboring cells
 
 ! calculation of direct neighbors
 ! x
 velgridcells(1) = act_cell
 IF(rel_pos(1) > 0) THEN
  velgridcells(2) = dyn_cell(act_cell)%neighbor(posx)
 ELSE
  velgridcells(2) = dyn_cell(act_cell)%neighbor(negx)
 END IF
 ! y
 IF(rel_pos(2) > 0) THEN
  velgridcells(3) = dyn_cell(act_cell)%neighbor(posy)
 ELSE
  velgridcells(3) = dyn_cell(act_cell)%neighbor(negy)
 END IF
 ! z
 IF(rel_pos(3) > 0) THEN
  velgridcells(4) = dyn_cell(act_cell)%neighbor(posz)
 ELSE
  velgridcells(4) = dyn_cell(act_cell)%neighbor(negz)
 END IF

 ! setting the testing vector pointing to the oposite cell in the block
 point_width = act_width
 point_corner = act_corner
 IF(rel_pos(1) < 0) THEN
  point_width(1) = - act_width(1)
  point_corner(1) = act_corner(1) + act_width(1)
 END IF
 IF(rel_pos(2) < 0) THEN
  point_width(2) = - act_width(2)
  point_corner(2) = act_corner(2) + act_width(2)
 END IF
 IF(rel_pos(3) < 0) THEN
  point_width(3) = - act_width(3)
  point_corner(3) = act_corner(3) + act_width(3)
 END IF
 point_vec = point_corner + 1.5*point_width
 write(28,*) act_corner, act_width
 write(28,*) point_corner, 1.5*point_width
 write(28,*) 0, 0, 0, point_vec

 CALL find_dyn_cell1(point_vec, point_cell)
 velgridcells(5) = point_cell

 ! neighbors of neighbors
 IF(rel_pos(1) > 0) THEN
  velgridcells(6) = dyn_cell(point_cell)%neighbor(negx)
 ELSE
  velgridcells(6) = dyn_cell(point_cell)%neighbor(posx)
 END IF
 ! y
 IF(rel_pos(2) > 0) THEN
  velgridcells(7) = dyn_cell(point_cell)%neighbor(negy)
 ELSE
  velgridcells(7) = dyn_cell(point_cell)%neighbor(posy)
 END IF
 ! z
 IF(rel_pos(3) > 0) THEN
  velgridcells(8) = dyn_cell(point_cell)%neighbor(negz)
 ELSE
  velgridcells(8) = dyn_cell(point_cell)%neighbor(posz)
 END IF

 DO I = 1,n_oct
  cur_cell = velgridcells(I)
  write(29,*) dyn_cell(cur_cell)%corner, dyn_cell(cur_cell)%width
 END DO
ELSE IF (dyngrid > 0) THEN
 STOP 'adaptive propagation grid is not currently suppotred'
! six possible  directions, distances to the boundaries and crosses
!  directions(1,:) = (/ 1, 0, 0 /)
!  distances(1) = act_corner(1) + act_width(1) - act_pos(1)
!  crossy(1) = posx
!  directions(2,:) = (/-1, 0, 0 /)
!  distances(2) = act_pos(1) - act_corner(1)
!  crossy(2) = negx
!  directions(3,:) = (/ 0, 1, 0 /)
!  distances(3) = act_corner(2) + act_width(2) - act_pos(2)
!  crossy(3) = posy
!  directions(4,:) = (/ 0,-1, 0 /)
!  distances(4) = act_pos(2) - act_corner(2)
!  crossy(4) = negy
!  directions(5,:) = (/ 0, 0, 1 /)
!  distances(5) = act_corner(3) + act_width(3) - act_pos(3)
!  crossy(5) = posz
!  directions(6,:) = (/ 0, 0,-1 /)
!  distances(6) = act_pos(3) - act_corner(3)
!  crossy(6) = negz
!  
!  DO cur_dir = 1,6
!   package(testPacket)%dir = directions(cur_dir, :)
!   package(testPacket)%next_cross = crossy(cur_dir)
!   CALL next_cell_down(testPacket, down_cell)
!   IF(dyngrid == 0) THEN
!    next_cell = down_cell
!   ELSE
!    CALL next_cell_up(testPacket, distances(cur_dir), down_cell, next_cell)
!   END IF
!   neighbors(cur_dir) = next_cell
!  END DO
END IF

! sort the velgridcells according to the index number
do J=2,n_oct
 I = J - 1

 a = velgridcells(J)

 do while(I >= 1)
  if(velgridcells(I) > a) THEN
   dummy = velgridcells(I+1)
   velgridcells(I+1) = velgridcells(I)
   velgridcells(I) = dummy
  end if
   I = I - 1
 end do
end do

END SUBROUTINE oct_neighbors
