! this subroutine will divide the given cell into eight smaller cells with
! half width
SUBROUTINE divide_cell_ijk(act_n_dyncell, max_n_dcell, numberofsubcells)

USE types

IMPLICIT NONE

! input variables
        ! number of created dynamical cells
        INTEGER                         :: max_n_dcell
        ! original dynamic cell
        INTEGER                         :: act_n_dyncell
        ! number of cells in the basic cells (nx, ny, nz)
        INTEGER, DIMENSION(3)           :: numberofsubcells
        ! properties of the original dynamic cell
        DOUBLE PRECISION, DIMENSION(3)  :: loc_corner, loc_cell_width
        DOUBLE PRECISION, DIMENSION(3)  :: new_cell_width
        ! loop indexes
        INTEGER                         :: I, J, K
        ! local cell index
        INTEGER                         :: L
        ! local number of subcells
        INTEGER                         :: loc_nx, loc_ny, loc_nz
        ! indexes of the next cells
        INTEGER                         :: xp, xm, yp, ym, zp, zm

dyn_cell(act_n_dyncell)%up_cell = max_n_dcell + 1
loc_corner = dyn_cell(act_n_dyncell)%corner
loc_cell_width = dyn_cell(act_n_dyncell)%width

loc_nx = numberofsubcells(1)
loc_ny = numberofsubcells(2)
loc_nz = numberofsubcells(3)

! definition of width of new cells
new_cell_width(1) = loc_cell_width(1) / DBLE(loc_nx)
new_cell_width(2) = loc_cell_width(2) / DBLE(loc_ny)
new_cell_width(3) = loc_cell_width(3) / DBLE(loc_nz)

print*, 'divide_cell_ijk: ', loc_nx, loc_ny, loc_nz

L = 1
DO I = 1, loc_nx
 DO J = 1, loc_ny
  DO K = 1, loc_nz
   dyn_cell(max_n_dcell + L)%width = loc_cell_width
   dyn_cell(max_n_dcell + L)%down_cell = act_n_dyncell
   dyn_cell(max_n_dcell + L)%up_cell = 0
   ! positions of newly created corners
   dyn_cell(max_n_dcell + L)%corner(1) = loc_corner(1) + DBLE((I - 1)) * new_cell_width(1)
   dyn_cell(max_n_dcell + L)%corner(2) = loc_corner(2) + DBLE((J - 1)) * new_cell_width(2)
   dyn_cell(max_n_dcell + L)%corner(3) = loc_corner(3) + DBLE((K - 1)) * new_cell_width(3)
   ! calculation of neighbors
   ! x+
   IF(I == nx_cell) THEN
    xp = 0
   ELSE
    xp = max_n_dcell + L + ny_cell * nz_cell
   END IF
   ! x-
   IF(I == 1) THEN
    xm = 0
   ELSE
    xm = max_n_dcell + L - ny_cell * nz_cell
   END IF
   ! y+
   IF(J == ny_cell) THEN
    yp = 0
   ELSE
    yp = max_n_dcell + L + nz_cell
   END IF
   ! y-
   IF(J == 1) THEN
    ym = 0
   ELSE
    ym = max_n_dcell + L - nz_cell
   END IF
   ! z+
   IF(K == nz_cell) THEN
    zp = 0
   ELSE
    zp = max_n_dcell + L + 1
   END IF
   ! z-
   IF(K == 1) THEN
    zm = 0
   ELSE
    zm = max_n_dcell + L - 1
   END IF
   ! setting number of neighbors
   dyn_cell(max_n_dcell + L)%neighbor(1) = xp
   dyn_cell(max_n_dcell + L)%neighbor(2) = xm
   dyn_cell(max_n_dcell + L)%neighbor(3) = yp
   dyn_cell(max_n_dcell + L)%neighbor(4) = ym
   dyn_cell(max_n_dcell + L)%neighbor(5) = zp
   dyn_cell(max_n_dcell + L)%neighbor(6) = zm
   L = L + 1
  END DO
 END DO
END DO

END SUBROUTINE divide_cell_ijk
