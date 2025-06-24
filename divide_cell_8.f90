! this subroutine will divide the given cell into eight smaller cells with
! half width
!
! INPUT: act_n_dyncell(INT): current index of a propGrid cell
! INTPUT/OUTPUT: max_n_dcell(INT): number of already created propGrid cells
!
SUBROUTINE divide_cell_8(act_n_dyncell, max_n_dcell)

USE types
USE constants

IMPLICIT NONE

! input variables
! number of created dynamical cells
INTEGER                         :: max_n_dcell
! original dynamic cell
INTEGER                         :: act_n_dyncell
INTEGER, PARAMETER              :: no_dcells = 8
INTEGER                         :: ind_I
! properties of the original dynamic cell
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: loc_corner, loc_cell_width
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: half_point, end_point

dyn_cell(act_n_dyncell)%up_cell = max_n_dcell + 1
loc_corner = dyn_cell(act_n_dyncell)%corner
loc_cell_width = dyn_cell(act_n_dyncell)%width
half_point = dyn_cell(act_n_dyncell)%corner + dyn_cell(act_n_dyncell)%width / 2.D0
end_point = dyn_cell(act_n_dyncell)%upcorner

DO ind_I = 1, no_dcells
 dyn_cell(max_n_dcell + ind_I)%width(ind_x) = (half_point(ind_x) - loc_corner(ind_x))
 dyn_cell(max_n_dcell + ind_I)%width(ind_y) = (half_point(ind_y) - loc_corner(ind_y))
 dyn_cell(max_n_dcell + ind_I)%width(ind_z) = (half_point(ind_z) - loc_corner(ind_z))
 dyn_cell(max_n_dcell + ind_I)%down_cell = act_n_dyncell
 dyn_cell(max_n_dcell + ind_I)%up_cell = 0
END DO
! the first cell
dyn_cell(max_n_dcell + 1)%corner(ind_x) = loc_corner(ind_x)
dyn_cell(max_n_dcell + 1)%corner(ind_y) = loc_corner(ind_y)
dyn_cell(max_n_dcell + 1)%corner(ind_z) = loc_corner(ind_z)

dyn_cell(max_n_dcell + 1)%upcorner(:) = half_point(:)
dyn_cell(max_n_dcell + 1)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 1)%neighbor(posx) = max_n_dcell + 2
dyn_cell(max_n_dcell + 1)%neighbor(negx) = 0
dyn_cell(max_n_dcell + 1)%neighbor(posy) = max_n_dcell + 3
dyn_cell(max_n_dcell + 1)%neighbor(negy) = 0
dyn_cell(max_n_dcell + 1)%neighbor(posz) = max_n_dcell + 5
dyn_cell(max_n_dcell + 1)%neighbor(negz) = 0
! the second cell
dyn_cell(max_n_dcell + 2)%corner(ind_x)  = half_point(ind_x)
dyn_cell(max_n_dcell + 2)%corner(ind_y)  = loc_corner(ind_y)
dyn_cell(max_n_dcell + 2)%corner(ind_z)  = loc_corner(ind_z)
dyn_cell(max_n_dcell + 2)%upcorner(ind_x)  = end_point(ind_x)
dyn_cell(max_n_dcell + 2)%upcorner(ind_y)  = loc_corner(ind_y)
dyn_cell(max_n_dcell + 2)%upcorner(ind_z)  = loc_corner(ind_z)
dyn_cell(max_n_dcell + 2)%down_cell      = act_n_dyncell
dyn_cell(max_n_dcell + 2)%neighbor(posx) = 0
dyn_cell(max_n_dcell + 2)%neighbor(negx) = max_n_dcell + 1
dyn_cell(max_n_dcell + 2)%neighbor(posy) = max_n_dcell + 4
dyn_cell(max_n_dcell + 2)%neighbor(negy) = 0
dyn_cell(max_n_dcell + 2)%neighbor(posz) = max_n_dcell + 6
dyn_cell(max_n_dcell + 2)%neighbor(negz) = 0
! the third cell
dyn_cell(max_n_dcell + 3)%corner(ind_x) = loc_corner(ind_x)
dyn_cell(max_n_dcell + 3)%corner(ind_y) = loc_corner(ind_y) + loc_cell_width(ind_y) / 2.E0
dyn_cell(max_n_dcell + 3)%corner(ind_z) = loc_corner(ind_z)
dyn_cell(max_n_dcell + 3)%upcorner(ind_x) = half_point(ind_x)
dyn_cell(max_n_dcell + 3)%upcorner(ind_y) = end_point(ind_y)
dyn_cell(max_n_dcell + 3)%upcorner(ind_z) = half_point(ind_z)
dyn_cell(max_n_dcell + 3)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 3)%neighbor(posx) = max_n_dcell + 4
dyn_cell(max_n_dcell + 3)%neighbor(negx) = 0
dyn_cell(max_n_dcell + 3)%neighbor(posy) = 0
dyn_cell(max_n_dcell + 3)%neighbor(negy) = max_n_dcell + 1
dyn_cell(max_n_dcell + 3)%neighbor(posz) = max_n_dcell + 7
dyn_cell(max_n_dcell + 3)%neighbor(negz) = 0
! the forth cell
dyn_cell(max_n_dcell + 4)%corner(ind_x) = loc_corner(ind_x) + loc_cell_width(ind_x) / 2.E0
dyn_cell(max_n_dcell + 4)%corner(ind_y) = loc_corner(ind_y) + loc_cell_width(ind_y) / 2.E0
dyn_cell(max_n_dcell + 4)%corner(ind_z) = loc_corner(ind_z)
dyn_cell(max_n_dcell + 4)%upcorner(ind_x) = end_point(ind_x)
dyn_cell(max_n_dcell + 4)%upcorner(ind_y) = end_point(ind_y)
dyn_cell(max_n_dcell + 4)%upcorner(ind_z) = half_point(ind_z)
dyn_cell(max_n_dcell + 4)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 4)%neighbor(posx) = 0
dyn_cell(max_n_dcell + 4)%neighbor(negx) = max_n_dcell + 3
dyn_cell(max_n_dcell + 4)%neighbor(posy) = 0
dyn_cell(max_n_dcell + 4)%neighbor(negy) = max_n_dcell + 2
dyn_cell(max_n_dcell + 4)%neighbor(posz) = max_n_dcell + 8
dyn_cell(max_n_dcell + 4)%neighbor(negz) = 0
! the fifth cell
dyn_cell(max_n_dcell + 5)%corner(ind_x) = loc_corner(ind_x)
dyn_cell(max_n_dcell + 5)%corner(ind_y) = loc_corner(ind_y)
dyn_cell(max_n_dcell + 5)%corner(ind_z) = loc_corner(ind_z) + loc_cell_width(ind_z) / 2.E0
dyn_cell(max_n_dcell + 5)%upcorner(ind_x) = half_point(ind_x)
dyn_cell(max_n_dcell + 5)%upcorner(ind_y) = half_point(ind_y)
dyn_cell(max_n_dcell + 5)%upcorner(ind_z) = end_point(ind_z)
dyn_cell(max_n_dcell + 5)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 5)%neighbor(posx) = max_n_dcell + 6
dyn_cell(max_n_dcell + 5)%neighbor(negx) = 0
dyn_cell(max_n_dcell + 5)%neighbor(posy) = max_n_dcell + 7
dyn_cell(max_n_dcell + 5)%neighbor(negy) = 0
dyn_cell(max_n_dcell + 5)%neighbor(posz) = 0
dyn_cell(max_n_dcell + 5)%neighbor(negz) = max_n_dcell + 1
! the sixth cell
dyn_cell(max_n_dcell + 6)%corner(ind_x) = loc_corner(ind_x) + loc_cell_width(ind_x) / 2.E0
dyn_cell(max_n_dcell + 6)%corner(ind_y) = loc_corner(ind_y)
dyn_cell(max_n_dcell + 6)%corner(ind_z) = loc_corner(ind_z) + loc_cell_width(ind_z) / 2.E0
dyn_cell(max_n_dcell + 6)%upcorner(ind_x) = end_point(ind_x)
dyn_cell(max_n_dcell + 6)%upcorner(ind_y) = half_point(ind_y)
dyn_cell(max_n_dcell + 6)%upcorner(ind_z) = end_point(ind_z)
dyn_cell(max_n_dcell + 6)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 6)%neighbor(posx) = 0
dyn_cell(max_n_dcell + 6)%neighbor(negx) = max_n_dcell + 5
dyn_cell(max_n_dcell + 6)%neighbor(posy) = max_n_dcell + 8
dyn_cell(max_n_dcell + 6)%neighbor(negy) = 0
dyn_cell(max_n_dcell + 6)%neighbor(posz) = 0
dyn_cell(max_n_dcell + 6)%neighbor(negz) = max_n_dcell + 2
! the seventh cell
dyn_cell(max_n_dcell + 7)%corner(ind_x) = loc_corner(ind_x)
dyn_cell(max_n_dcell + 7)%corner(ind_y) = loc_corner(ind_y) + loc_cell_width(ind_y) / 2.E0
dyn_cell(max_n_dcell + 7)%corner(ind_z) = loc_corner(ind_z) + loc_cell_width(ind_z) / 2.E0
dyn_cell(max_n_dcell + 7)%upcorner(ind_x) = half_point(ind_x)
dyn_cell(max_n_dcell + 7)%upcorner(ind_y) = end_point(ind_y)
dyn_cell(max_n_dcell + 7)%upcorner(ind_z) = end_point(ind_z)
dyn_cell(max_n_dcell + 7)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 7)%neighbor(posx) = max_n_dcell + 8
dyn_cell(max_n_dcell + 7)%neighbor(negx) = 0
dyn_cell(max_n_dcell + 7)%neighbor(posy) = 0
dyn_cell(max_n_dcell + 7)%neighbor(negy) = max_n_dcell + 5
dyn_cell(max_n_dcell + 7)%neighbor(posz) = 0
dyn_cell(max_n_dcell + 7)%neighbor(negz) = max_n_dcell + 3
! the eighth cell
dyn_cell(max_n_dcell + no_dcells)%corner(ind_x) = loc_corner(ind_x) + loc_cell_width(ind_x) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%corner(ind_y) = loc_corner(ind_y) + loc_cell_width(ind_y) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%corner(ind_z) = loc_corner(ind_z) + loc_cell_width(ind_z) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%upcorner(ind_x) = end_point(ind_x)
dyn_cell(max_n_dcell + no_dcells)%upcorner(ind_y) = end_point(ind_y)
dyn_cell(max_n_dcell + no_dcells)%upcorner(ind_z) = end_point(ind_z)
dyn_cell(max_n_dcell + no_dcells)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + no_dcells)%neighbor(posx) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbor(negx) = max_n_dcell + 7
dyn_cell(max_n_dcell + no_dcells)%neighbor(posy) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbor(negy) = max_n_dcell + 6
dyn_cell(max_n_dcell + no_dcells)%neighbor(posz) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbor(negz) = max_n_dcell + 4

END SUBROUTINE divide_cell_8
