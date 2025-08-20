! finds eight points for trilinear interpolation
!
! INPUT: pack_index(INT) -- index of package
!        rel_pos(DBLE(const_dimofspace)) -- relative position of packet in respect to the propGrid cell centre
! OUTPUT: cube_pos(DBLE(8, const_dimofspace)) -- position of interpolation points
!         incell(LOG) -- is the current propGrid cell index < 0
!
SUBROUTINE oct_virtcube(pack_index, cube_pos, velgridcells, incell)

USE types
USE constants
USE dummypacket
IMPLICIT NONE

INTEGER                                                 :: pack_index, dummypack_index

INTEGER, PARAMETER                                      :: n_oct=8

DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: rel_pos
LOGICAL, DIMENSION(const_dimofspace)                                   :: posxyz
INTEGER                                                 :: sgn_x, sgn_y, sgn_z
DOUBLE PRECISION, DIMENSION(8,const_dimofspace)                        :: cube_pos

INTEGER                                                 :: act_cell, cur_cell
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: act_pos, act_corner, act_width, act_center, act_upcorner
DOUBLE PRECISION, DIMENSION(const_dimofspace)                          :: cur_pos
INTEGER, DIMENSION(n_oct)                   :: velgridcells
INTEGER, PARAMETER                      :: n_zero = 1, n_x = 2, n_y = 3, n_z = 4,&
                                           n_xz = 5, n_xy = 6, n_xyz = 7, n_yz = 8

LOGICAL                                                 :: incell
INTEGER, PARAMETER                      :: dir_x = 1, dir_y = 2, dir_z = 3
INTEGER                                                 :: ind_I

IF(pack_index <= SIZE(package)) THEN
 act_cell = package(pack_index)%cell_numb
 act_pos = package(pack_index)%pos
ELSE IF(pack_index > SIZE(package)) THEN
 dummypack_index = pack_index - SIZE(package)
 act_cell = dummypackage(dummypack_index)%cell_numb
 act_pos = dummypackage(dummypack_index)%pos
END IF
act_corner = dyn_cell(act_cell)%corner
act_upcorner = dyn_cell(act_cell)%upcorner
act_width = act_upcorner - act_corner
act_center = (act_corner + act_upcorner)/2.0

rel_pos = act_pos - act_center

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
DO ind_I = 1,8
 cube_pos(ind_I,:) = act_center(:)
END DO

! z
cube_pos(2, ind_z) = act_pos(ind_z) + sgn_z * act_width(ind_z)
! y
cube_pos(3, ind_y) = act_pos(ind_y) + sgn_y * act_width(ind_y)
! zy
cube_pos(4, ind_y) = act_pos(ind_y) + sgn_y * act_width(ind_y)
cube_pos(4, ind_z) = act_pos(ind_z) + sgn_z * act_width(ind_z)
! x
cube_pos(5, ind_x) = act_pos(ind_x) + sgn_x * act_width(ind_x)
! zx
cube_pos(6, ind_x) = act_pos(ind_x) + sgn_x * act_width(ind_x)
cube_pos(6, ind_z) = act_pos(ind_z) + sgn_z * act_width(ind_z)
! xy
cube_pos(7, ind_x) = act_pos(ind_x) + sgn_x * act_width(ind_x)
cube_pos(7, ind_y) = act_pos(ind_y) + sgn_y * act_width(ind_y)
! zxy
cube_pos(8, ind_x) = act_pos(ind_x) + sgn_x * act_width(ind_x)
cube_pos(8, ind_y) = act_pos(ind_y) + sgn_y * act_width(ind_y)
cube_pos(8, ind_z) = act_pos(ind_z) + sgn_z * act_width(ind_z)

velgridcells(ind_x) = act_cell
! the current propagation grid cells
DO ind_I = 2, 8
 cur_pos = cube_pos(ind_I,:)
 CALL find_dyn_cell1(cur_pos, cur_cell)
 velgridcells(ind_I) = cur_cell
 IF(cur_cell < 0) THEN
  incell = .TRUE.
  RETURN
 END IF
END DO

! DO ind_I = 1,8
!  write(29,*) cube_pos(I,:)
! END DO


! sort the velgridcells according to the index number
! do J = 2, n_oct
!  ind_I = J - 1
! 
!  a = velgridcells(J)
! 
!  do while(ind_I >= 1)
!   if(velgridcells(ind_I) > a) THEN
!    dummy = velgridcells(ind_I + 1)
!    dummy2 = cube_pos(ind_I + 1, :)
!    velgridcells(ind_I + 1) = velgridcells(I)
!    cube_pos(ind_I + 1, :) = cube_pos(I, :)
!    velgridcells(ind_I) = dummy
!    cube_pos(ind_I, :) = dummy2
!   end if
!    ind_I = ind_I - 1
!  end do
! end do

! STOP 'oct_virtcube: testing'

END SUBROUTINE oct_virtcube
