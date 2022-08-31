SUBROUTINE oct_neighbors(pack_index, rel_pos, velgridcells, incell)


USE types
IMPLICIT NONE

INTEGER, PARAMETER                      :: n_oct=8

INTEGER                                 :: pack_index
DOUBLE PRECISION, DIMENSION(3)          :: rel_pos

INTEGER                                 :: act_cell
DOUBLE PRECISION, DIMENSION(3)          :: act_corner, act_width
DOUBLE PRECISION, DIMENSION(3)          :: act_center

INTEGER                                 :: cur_cell

INTEGER, DIMENSION(n_oct)                   :: velgridcells

INTEGER                                 :: I, J, dummy, a

LOGICAL                                 :: incell
LOGICAL, DIMENSION(3)                   :: posxyz

INTEGER, PARAMETER                      :: n_zero = 1, n_x = 2, n_y = 3, n_z = 4,&
                                           n_xz = 5, n_xy = 6, n_xyz = 7, n_yz = 8
INTEGER, PARAMETER                      :: dir_x = 1, dir_y = 2, dir_z = 3
INTEGER                                 :: neigx, neigy
INTEGER                                 :: neigxy, neigxz, neigxyz, neigyz

DOUBLE PRECISION, DIMENSION(3)          :: act_pos


act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width
act_center = act_corner + act_width/2.0

! write(*,*) 'oct_neighbors: rel_pos = ', rel_pos
! 
! write(24,*) act_corner, act_width

! this part code only for the regular grid
! posxyz(1-3) indicates a move in the three directions from the zero point
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


! test for the boundary propagation cell
DO I = 1,n_oct
 cur_cell = velgridcells(I)
 incell = .FALSE.
 IF(cur_cell < 0) THEN
  incell = .TRUE.
  EXIT
 END IF
END DO
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
