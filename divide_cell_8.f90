! this subroutine will divide the given cell into eight smaller cells with
! half width
SUBROUTINE divide_cell_8(act_n_dyncell, max_n_dcell)

USE types

IMPLICIT NONE

! input variables
        ! number of created dynamical cells
        INTEGER                         :: max_n_dcell
        ! original dynamic cell
        INTEGER                         :: act_n_dyncell
        INTEGER, PARAMETER              :: no_dcells = 8
        INTEGER                         :: I
        ! properties of the original dynamic cell
        DOUBLE PRECISION, DIMENSION(3)  :: loc_corner, loc_cell_width


loc_corner(1) = dyn_cell(act_n_dyncell)%corner(1)
loc_corner(2) = dyn_cell(act_n_dyncell)%corner(2)
loc_corner(3) = dyn_cell(act_n_dyncell)%corner(3)
loc_cell_width(1) = dyn_cell(act_n_dyncell)%width(1)
loc_cell_width(2) = dyn_cell(act_n_dyncell)%width(2)
loc_cell_width(3) = dyn_cell(act_n_dyncell)%width(3)

DO I = 1, no_dcells
 dyn_cell(max_n_dcell + I)%width(1) = loc_cell_width(1) / 2.D0
 dyn_cell(max_n_dcell + I)%width(2) = loc_cell_width(2) / 2.D0
 dyn_cell(max_n_dcell + I)%width(3) = loc_cell_width(3) / 2.D0
 dyn_cell(max_n_dcell + I)%down_cell = act_n_dyncell
 dyn_cell(max_n_dcell + I)%up_cell = 0
END DO
! the first cell
dyn_cell(max_n_dcell + 1)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 1)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 1)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 1)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 1)%neighbor(1) = max_n_dcell + 2
dyn_cell(max_n_dcell + 1)%neighbor(2) = 0
dyn_cell(max_n_dcell + 1)%neighbor(3) = max_n_dcell + 3
dyn_cell(max_n_dcell + 1)%neighbor(4) = 0
dyn_cell(max_n_dcell + 1)%neighbor(5) = max_n_dcell + 5
dyn_cell(max_n_dcell + 1)%neighbor(6) = 0
! the second cell
dyn_cell(max_n_dcell + 2)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + 2)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 2)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 2)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 2)%neighbor(1) = 0
dyn_cell(max_n_dcell + 2)%neighbor(2) = max_n_dcell + 1
dyn_cell(max_n_dcell + 2)%neighbor(3) = max_n_dcell + 4
dyn_cell(max_n_dcell + 2)%neighbor(4) = 0
dyn_cell(max_n_dcell + 2)%neighbor(5) = max_n_dcell + 6
dyn_cell(max_n_dcell + 2)%neighbor(6) = 0
! the third cell
dyn_cell(max_n_dcell + 3)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 3)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + 3)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 3)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 3)%neighbor(1) = max_n_dcell + 4
dyn_cell(max_n_dcell + 3)%neighbor(2) = 0
dyn_cell(max_n_dcell + 3)%neighbor(3) = 0
dyn_cell(max_n_dcell + 3)%neighbor(4) = max_n_dcell + 1
dyn_cell(max_n_dcell + 3)%neighbor(5) = max_n_dcell + 7
dyn_cell(max_n_dcell + 3)%neighbor(6) = 0
! the forth cell
dyn_cell(max_n_dcell + 4)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + 4)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + 4)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 4)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 4)%neighbor(1) = 0
dyn_cell(max_n_dcell + 4)%neighbor(2) = max_n_dcell + 3
dyn_cell(max_n_dcell + 4)%neighbor(3) = 0
dyn_cell(max_n_dcell + 4)%neighbor(4) = max_n_dcell + 2
dyn_cell(max_n_dcell + 4)%neighbor(5) = max_n_dcell + 8
dyn_cell(max_n_dcell + 4)%neighbor(6) = 0
! the fifth cell
dyn_cell(max_n_dcell + 5)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 5)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 5)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + 5)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 5)%neighbor(1) = max_n_dcell + 6
dyn_cell(max_n_dcell + 5)%neighbor(2) = 0
dyn_cell(max_n_dcell + 5)%neighbor(3) = max_n_dcell + 7
dyn_cell(max_n_dcell + 5)%neighbor(4) = 0
dyn_cell(max_n_dcell + 5)%neighbor(5) = 0
dyn_cell(max_n_dcell + 5)%neighbor(6) = max_n_dcell + 1
! the sixth cell
dyn_cell(max_n_dcell + 6)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + 6)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 6)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + 6)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 6)%neighbor(1) = 0
dyn_cell(max_n_dcell + 6)%neighbor(2) = max_n_dcell + 5
dyn_cell(max_n_dcell + 6)%neighbor(3) = max_n_dcell + 8
dyn_cell(max_n_dcell + 6)%neighbor(4) = 0
dyn_cell(max_n_dcell + 6)%neighbor(5) = 0
dyn_cell(max_n_dcell + 6)%neighbor(6) = max_n_dcell + 2
! the seventh cell
dyn_cell(max_n_dcell + 7)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 7)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + 7)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + 7)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 7)%neighbor(1) = max_n_dcell + 8
dyn_cell(max_n_dcell + 7)%neighbor(2) = 0
dyn_cell(max_n_dcell + 7)%neighbor(3) = 0
dyn_cell(max_n_dcell + 7)%neighbor(4) = max_n_dcell + 5
dyn_cell(max_n_dcell + 7)%neighbor(5) = 0
dyn_cell(max_n_dcell + 7)%neighbor(6) = max_n_dcell + 3
! the eighth cell
dyn_cell(max_n_dcell + no_dcells)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + no_dcells)%neighbor(1) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbor(2) = max_n_dcell + 7
dyn_cell(max_n_dcell + no_dcells)%neighbor(3) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbor(4) = max_n_dcell + 6
dyn_cell(max_n_dcell + no_dcells)%neighbor(5) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbor(6) = max_n_dcell + 4

END SUBROUTINE divide_cell_8
