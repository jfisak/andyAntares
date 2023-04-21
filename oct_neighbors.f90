SUBROUTINE oct_neighbors(pack_index, rel_pos, velgridcells, incell)


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

LOGICAL                                 :: incell
LOGICAL, DIMENSION(3)                   :: posxyz


INTEGER, PARAMETER                      :: n_zero = 1, n_x = 2, n_y = 3, n_z = 4,&
                                           n_xz = 5, n_xy = 6, n_xyz = 7, n_yz = 8
INTEGER, PARAMETER                      :: dir_x = 1, dir_y = 2, dir_z = 3
INTEGER                                 :: neigx, neigy, neigz
INTEGER                                 :: neigxy, neigxz, neigxyz, neigyz

act_cell = package(pack_index)%cell_numb
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width

! write(*,*) 'oct_neighbors: rel_pos = ', rel_pos
! 
! write(24,*) act_corner, act_width

! this part code only for the regular grid
! posxyz(1-3) indicates a move in the three directions from the zero point
IF(dyngrid == 0) THEN
 ! find the neighboring cells
 
 ! calculation of direct neighbors
 ! x
 velgridcells(n_zero) = act_cell
 IF(rel_pos(dir_x) > 0) THEN
  velgridcells(n_x) = dyn_cell(act_cell)%neighbor(posx)
  posxyz(dir_x) = .TRUE.
 ELSE
  velgridcells(n_x) = dyn_cell(act_cell)%neighbor(negx)
  posxyz(dir_x) = .FALSE.
 END IF
 ! y
 IF(rel_pos(dir_y) > 0) THEN
  velgridcells(n_y) = dyn_cell(act_cell)%neighbor(posy)
  posxyz(dir_y) = .TRUE.
 ELSE
  velgridcells(n_y) = dyn_cell(act_cell)%neighbor(negy)
  posxyz(dir_y) = .FALSE.
 END IF
 ! z
 IF(rel_pos(dir_z) > 0) THEN
  velgridcells(n_z) = dyn_cell(act_cell)%neighbor(posz)
  posxyz(dir_z) = .TRUE.
 ELSE
  velgridcells(n_z) = dyn_cell(act_cell)%neighbor(negz)
  posxyz(dir_z) = .FALSE.
 END IF
 ! test if the cell is on the edge of the propGrid
 DO I = 2,4
  IF(velgridcells(I) < 0) THEN
   incell = .true.
   RETURN
  END IF
 END DO
 
 ! neighbors x
 neigx = velgridcells(n_x)
 ! neighbor xz
 ! neighbor xy
 IF(posxyz(dir_y)) THEN
  neigxy = dyn_cell(neigx)%neighbor(posy)
 ELSE
  neigxy = dyn_cell(neigx)%neighbor(negy)
 END IF
 velgridcells(n_xy) = neigxy
 ! neighbor xz
 ! neighbor xyz
 IF(posxyz(dir_z)) THEN
  neigxz = dyn_cell(neigx)%neighbor(posz)
  neigxyz = dyn_cell(neigxy)%neighbor(posz)
 ELSE
  neigxz = dyn_cell(neigx)%neighbor(negz)
  neigxyz = dyn_cell(neigxy)%neighbor(negz)
 END IF
 velgridcells(n_xz)=neigxz
 velgridcells(n_xyz)=neigxyz

 ! neighbors y
 neigy = velgridcells(n_y)
 ! neighbor yz
 IF(posxyz(dir_z)) THEN
  neigyz = dyn_cell(neigy)%neighbor(posz)
 ELSE
  neigyz = dyn_cell(neigy)%neighbor(negz)
 END IF
 velgridcells(n_yz) = neigyz

 ! write(*,*) 'oct_neighbors: velgridcells = ', velgridcells
 

 DO I = 1,n_oct
  cur_cell = velgridcells(I)
  ! IF(cur_cell > 0) write(29,*) dyn_cell(cur_cell)%corner, dyn_cell(cur_cell)%width
  incell = .FALSE.
  IF(cur_cell < 0) THEN
   incell = .TRUE.
   EXIT
  END IF
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

! DO I = 1,8
!  IF(velgridcells(I) < 0) STOP 'oct_neighbors: index < 0'
! END DO


END SUBROUTINE oct_neighbors
