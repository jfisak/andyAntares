SUBROUTINE create_subcells(max_n_dcell, act_n_dyncell, no_dcells)

USE types

IMPLICIT NONE
   
! input varibles
INTEGER                                 :: max_n_dcell, act_n_dyncell, no_dcells
DOUBLE PRECISION, DIMENSION(3)          :: loc_corner, loc_cell_width
INTEGER                                 :: I

loc_corner(1) = dyn_cell(act_n_dyncell)%corner(1)
loc_corner(2) = dyn_cell(act_n_dyncell)%corner(2)
loc_corner(3) = dyn_cell(act_n_dyncell)%corner(3)
loc_cell_width(1) = dyn_cell(act_n_dyncell)%width(1)
loc_cell_width(2) = dyn_cell(act_n_dyncell)%width(3)
loc_cell_width(3) = dyn_cell(act_n_dyncell)%width(2)

DO I = 1, no_dcells
 dyn_cell(max_n_dcell + I)%width(1) = loc_cell_width(1) / 2.E0
 dyn_cell(max_n_dcell + I)%width(2) = loc_cell_width(2) / 2.E0
 dyn_cell(max_n_dcell + I)%width(3) = loc_cell_width(3) / 2.E0
 dyn_cell(max_n_dcell + I)%down_cell = act_n_dyncell
 dyn_cell(max_n_dcell + I)%up_cell = 0
END DO
! the first cell
dyn_cell(max_n_dcell + 1)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 1)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 1)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 1)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 1)%neighbour(1) = max_n_dcell + 2
dyn_cell(max_n_dcell + 1)%neighbour(2) = 0
dyn_cell(max_n_dcell + 1)%neighbour(3) = max_n_dcell + 3
dyn_cell(max_n_dcell + 1)%neighbour(4) = 0
dyn_cell(max_n_dcell + 1)%neighbour(5) = max_n_dcell + 5
dyn_cell(max_n_dcell + 1)%neighbour(6) = 0
! the second cell
dyn_cell(max_n_dcell + 2)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + 2)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 2)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 2)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 2)%neighbour(1) = 0
dyn_cell(max_n_dcell + 2)%neighbour(2) = max_n_dcell + 1
dyn_cell(max_n_dcell + 2)%neighbour(3) = max_n_dcell + 4
dyn_cell(max_n_dcell + 2)%neighbour(4) = 0
dyn_cell(max_n_dcell + 2)%neighbour(5) = max_n_dcell + 6
dyn_cell(max_n_dcell + 2)%neighbour(6) = 0
! the third cell
dyn_cell(max_n_dcell + 3)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 3)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + 3)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 3)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 3)%neighbour(1) = max_n_dcell + 4
dyn_cell(max_n_dcell + 3)%neighbour(2) = 0
dyn_cell(max_n_dcell + 3)%neighbour(3) = 0
dyn_cell(max_n_dcell + 3)%neighbour(4) = max_n_dcell + 1
dyn_cell(max_n_dcell + 3)%neighbour(5) = max_n_dcell + 7
dyn_cell(max_n_dcell + 3)%neighbour(6) = 0
! the forth cell
dyn_cell(max_n_dcell + 4)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + 4)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + 4)%corner(3) = loc_corner(3)
dyn_cell(max_n_dcell + 4)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 4)%neighbour(1) = 0
dyn_cell(max_n_dcell + 4)%neighbour(2) = max_n_dcell + 3
dyn_cell(max_n_dcell + 4)%neighbour(3) = 0
dyn_cell(max_n_dcell + 4)%neighbour(4) = max_n_dcell + 2
dyn_cell(max_n_dcell + 4)%neighbour(5) = max_n_dcell + 8
dyn_cell(max_n_dcell + 4)%neighbour(6) = 0
! the fifth cell
dyn_cell(max_n_dcell + 5)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 5)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 5)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + 5)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 5)%neighbour(1) = max_n_dcell + 6
dyn_cell(max_n_dcell + 5)%neighbour(2) = 0
dyn_cell(max_n_dcell + 5)%neighbour(3) = max_n_dcell + 7
dyn_cell(max_n_dcell + 5)%neighbour(4) = 0
dyn_cell(max_n_dcell + 5)%neighbour(5) = 0
dyn_cell(max_n_dcell + 5)%neighbour(6) = max_n_dcell + 1
! the sixth cell
dyn_cell(max_n_dcell + 6)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + 6)%corner(2) = loc_corner(2)
dyn_cell(max_n_dcell + 6)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + 6)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 6)%neighbour(1) = 0
dyn_cell(max_n_dcell + 6)%neighbour(2) = max_n_dcell + 5
dyn_cell(max_n_dcell + 6)%neighbour(3) = max_n_dcell + 8
dyn_cell(max_n_dcell + 6)%neighbour(4) = 0
dyn_cell(max_n_dcell + 6)%neighbour(5) = 0
dyn_cell(max_n_dcell + 6)%neighbour(6) = max_n_dcell + 2
! the seventh cell
dyn_cell(max_n_dcell + 7)%corner(1) = loc_corner(1)
dyn_cell(max_n_dcell + 7)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + 7)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + 7)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + 7)%neighbour(1) = max_n_dcell + 8
dyn_cell(max_n_dcell + 7)%neighbour(2) = 0
dyn_cell(max_n_dcell + 7)%neighbour(3) = 0
dyn_cell(max_n_dcell + 7)%neighbour(4) = max_n_dcell + 5
dyn_cell(max_n_dcell + 7)%neighbour(5) = 0
dyn_cell(max_n_dcell + 7)%neighbour(6) = max_n_dcell + 3
! the eighth cell
dyn_cell(max_n_dcell + no_dcells)%corner(1) = loc_corner(1) + loc_cell_width(1) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%corner(2) = loc_corner(2) + loc_cell_width(2) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%corner(3) = loc_corner(3) + loc_cell_width(3) / 2.E0
dyn_cell(max_n_dcell + no_dcells)%down_cell = act_n_dyncell
dyn_cell(max_n_dcell + no_dcells)%neighbour(1) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbour(2) = max_n_dcell + 7
dyn_cell(max_n_dcell + no_dcells)%neighbour(3) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbour(4) = max_n_dcell + 6
dyn_cell(max_n_dcell + no_dcells)%neighbour(5) = 0
dyn_cell(max_n_dcell + no_dcells)%neighbour(6) = max_n_dcell + 4

END SUBROUTINE create_subcells
