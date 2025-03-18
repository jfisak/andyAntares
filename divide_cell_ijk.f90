! this subroutine will divide the given cell into a subgrid nx * ny * nz
!
! INPUT: act_n_dyncell(INT): current propGrid cell index
!        max_n_dcell(INT): maximun number of propGrid cells
!        numberofsubcells(INT): number of subcells
! OUTPUT: NONE
!
SUBROUTINE divide_cell_ijk(act_n_dyncell, max_n_dcell, numberofsubcells)

USE types
USE constants

IMPLICIT NONE

! input variables
! number of created dynamical cells
INTEGER                         :: max_n_dcell
! original dynamic cell
INTEGER                         :: act_n_dyncell
! number of cells in the basic cells (nx, ny, nz)
INTEGER, DIMENSION(const_dimofspace)           :: numberofsubcells
! properties of the original dynamic cell
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: loc_corner, loc_cell_width
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: new_cell_width
! loop indexes
INTEGER                         :: ind_I, ind_J, ind_K
! local cell index
INTEGER                         :: ind_L
! local number of subcells
INTEGER                         :: loc_nx, loc_ny, loc_nz
! indexes of the next cells
INTEGER                         :: xp, xm, yp, ym, zp, zm

dyn_cell(act_n_dyncell)%up_cell = max_n_dcell + 1
loc_corner = dyn_cell(act_n_dyncell)%corner
loc_cell_width = dyn_cell(act_n_dyncell)%width

loc_nx = numberofsubcells(ind_x)
loc_ny = numberofsubcells(ind_y)
loc_nz = numberofsubcells(ind_z)

! definition of width of new cells
new_cell_width(ind_x) = loc_cell_width(ind_x) / DBLE(loc_nx)
new_cell_width(ind_y) = loc_cell_width(ind_y) / DBLE(loc_ny)
new_cell_width(ind_z) = loc_cell_width(ind_z) / DBLE(loc_nz)

! write(*,*) 'divide_cell_ijk: r+n*w = ', loc_corner(1) + dble(loc_nx) * new_cell_width(1)
! write(*,*) 'divide_cell_ijk: BG r+w = ', dyn_cell(act_n_dyncell)%corner + dyn_cell(act_n_dyncell)%width

ind_L = 1
DO ind_I = 1, loc_nx
 DO ind_J = 1, loc_ny
  DO ind_K = 1, loc_nz
   dyn_cell(max_n_dcell + ind_L)%width = new_cell_width
   dyn_cell(max_n_dcell + ind_L)%down_cell = act_n_dyncell
   dyn_cell(max_n_dcell + ind_L)%up_cell = 0
   ! positions of newly created corners
   dyn_cell(max_n_dcell + ind_L)%corner(ind_x) = loc_corner(ind_x) + DBLE((ind_I - 1)) * new_cell_width(ind_x)
   dyn_cell(max_n_dcell + ind_L)%corner(ind_y) = loc_corner(ind_y) + DBLE((ind_J - 1)) * new_cell_width(ind_y)
   dyn_cell(max_n_dcell + ind_L)%corner(ind_z) = loc_corner(ind_z) + DBLE((ind_K - 1)) * new_cell_width(ind_z)
   IF(ind_I == loc_nx) THEN
    dyn_cell(max_n_dcell + ind_L)%width(ind_x) = dyn_cell(act_n_dyncell)%corner(ind_x) &
     + dyn_cell(act_n_dyncell)%width(ind_x) - dyn_cell(max_n_dcell + ind_L)%corner(ind_x)
   ELSE IF (ind_J == loc_ny) THEN
    dyn_cell(max_n_dcell + ind_L)%width(ind_y) = dyn_cell(act_n_dyncell)%corner(ind_y) &
     + dyn_cell(act_n_dyncell)%width(ind_y) - dyn_cell(max_n_dcell + ind_L)%corner(ind_y)
   ELSE IF (ind_K == loc_nz) THEN
    dyn_cell(max_n_dcell + ind_L)%width(ind_z) = dyn_cell(act_n_dyncell)%corner(ind_z) &
     + dyn_cell(act_n_dyncell)%width(ind_z) - dyn_cell(max_n_dcell + ind_L)%corner(ind_z)
   END IF
   ! calculation of neighbors
   ! x+
   IF(ind_I == loc_nx) THEN
    xp = 0
   ELSE
    xp = max_n_dcell + ind_L + loc_ny * loc_nz
   END IF
   ! x-
   IF(ind_I == 1) THEN
    xm = 0
   ELSE
    xm = max_n_dcell + ind_L - loc_ny * loc_nz
   END IF
   ! y+
   IF(ind_J == loc_ny) THEN
    yp = 0
   ELSE
    yp = max_n_dcell + ind_L + loc_nz
   END IF
   ! y-
   IF(ind_J == 1) THEN
    ym = 0
   ELSE
    ym = max_n_dcell + ind_L - loc_nz
   END IF
   ! z+
   IF(ind_K == loc_nz) THEN
    zp = 0
   ELSE
    zp = max_n_dcell + ind_L + 1
   END IF
   ! z-
   IF(ind_K == 1) THEN
    zm = 0
   ELSE
    zm = max_n_dcell + ind_L - 1
   END IF
   ! setting number of neighbors
   dyn_cell(max_n_dcell + ind_L)%neighbor(posx) = xp
   dyn_cell(max_n_dcell + ind_L)%neighbor(negx) = xm
   dyn_cell(max_n_dcell + ind_L)%neighbor(posy) = yp
   dyn_cell(max_n_dcell + ind_L)%neighbor(negy) = ym
   dyn_cell(max_n_dcell + ind_L)%neighbor(posz) = zp
   dyn_cell(max_n_dcell + ind_L)%neighbor(negz) = zm
   ind_L = ind_L + 1
   ! IF(ind_K == loc_nz) THEN
   !  write(*,*) 'divide_cell_ijk: K = ', K
   !  write(*,*) 'divide_cell_ijk: r + w = ', dyn_cell(max_n_dcell + ind_L - 1)%corner + dyn_cell(max_n_dcell + ind_L)%width
   !  write(*,*) 'divide_cell_ijk BC: r + w = ', dyn_cell(act_n_dyncell)%corner + dyn_cell(act_n_dyncell)%width
   ! END IF
  END DO
 END DO
END DO

END SUBROUTINE divide_cell_ijk
