SUBROUTINE oct_virtcube(pack_index, rel_pos, cube_pos, incell)

USE types
USE constants
IMPLICIT NONE

INTEGER                                                 :: pack_index

INTEGER, PARAMETER                                      :: n_oct=8

DOUBLE PRECISION, DIMENSION(3)                          :: rel_pos
LOGICAL, DIMENSION(3)                                   :: posxyz
INTEGER                                                 :: sgn_x, sgn_y, sgn_z
DOUBLE PRECISION, DIMENSION(8,3)                        :: cube_pos

INTEGER                                                 :: act_cell, cur_cell
DOUBLE PRECISION, DIMENSION(3)                          :: act_pos, act_corner, act_width, act_center
DOUBLE PRECISION, DIMENSION(3)                          :: cur_pos
INTEGER, DIMENSION(n_oct)                   :: velgridcells
INTEGER, PARAMETER                      :: n_zero = 1, n_x = 2, n_y = 3, n_z = 4,&
                                           n_xz = 5, n_xy = 6, n_xyz = 7, n_yz = 8

LOGICAL                                                 :: incell
INTEGER, PARAMETER                      :: dir_x = 1, dir_y = 2, dir_z = 3
INTEGER                                 :: I, J, dummy, a
DOUBLE PRECISION, DIMENSION(3)          :: dummy2

act_cell = package(pack_index)%cell_numb
act_pos = package(pack_index)%pos
act_corner = dyn_cell(act_cell)%corner
act_width = dyn_cell(act_cell)%width
act_center = act_corner + act_width/2.0

! x
IF(rel_pos(dir_x) > 0) THEN
 posxyz(dir_x) = .TRUE.
 sgn_x = 1
ELSE
 posxyz(dir_x) = .FALSE.
 sgn_x = -1
END IF

! y
IF(rel_pos(dir_y) > 0) THEN
 posxyz(dir_y) = .TRUE.
 sgn_y = 1
ELSE
 posxyz(dir_y) = .FALSE.
 sgn_y = -1
END IF

! z
IF(rel_pos(dir_z) > 0) THEN
 posxyz(dir_z) = .TRUE.
 sgn_z = 1
ELSE
 posxyz(dir_z) = .FALSE.
 sgn_z = -1
END IF

! we will create a regular cube with one point located in the
! center of the current cell and the other points in a current
! quadrant
! first of all positions of the cube
DO I = 1,8
 cube_pos(I,:) = act_center(:)
END DO

! z
cube_pos(2,3) = act_center(3) + sgn_z * act_width(3)
! y
cube_pos(3,2) = act_center(2) + sgn_y * act_width(2)
! zy
cube_pos(4,2) = act_center(2) + sgn_y * act_width(2)
cube_pos(4,3) = act_center(3) + sgn_z * act_width(3)
! x
cube_pos(5,1) = act_center(1) + sgn_x * act_width(1)
! zx
cube_pos(6,1) = act_center(1) + sgn_x * act_width(1)
cube_pos(6,3) = act_center(3) + sgn_z * act_width(3)
! xy
cube_pos(7,1) = act_center(1) + sgn_x * act_width(1)
cube_pos(7,2) = act_center(2) + sgn_y * act_width(2)
! zxy
cube_pos(8,1) = act_center(1) + sgn_x * act_width(1)
cube_pos(8,2) = act_center(2) + sgn_y * act_width(2)
cube_pos(8,3) = act_center(3) + sgn_z * act_width(3)

velgridcells(1) = act_cell
! the current propagation grid cells
DO I = 2, 8
 cur_pos = cube_pos(I,:)
 CALL find_dyn_cell1(cur_pos, cur_cell)
 velgridcells(I) = cur_cell
 IF(cur_cell < 0) THEN
  incell = .TRUE.
  RETURN
 END IF
END DO

! DO I = 1,8
!  write(29,*) cube_pos(I,:)
! END DO


! sort the velgridcells according to the index number
! do J = 2, n_oct
!  I = J - 1
! 
!  a = velgridcells(J)
! 
!  do while(I >= 1)
!   if(velgridcells(I) > a) THEN
!    dummy = velgridcells(I + 1)
!    dummy2 = cube_pos(I + 1, :)
!    velgridcells(I + 1) = velgridcells(I)
!    cube_pos(I + 1, :) = cube_pos(I, :)
!    velgridcells(I) = dummy
!    cube_pos(I, :) = dummy2
!   end if
!    I = I - 1
!  end do
! end do

! STOP 'oct_virtcube: testing'

END SUBROUTINE oct_virtcube
