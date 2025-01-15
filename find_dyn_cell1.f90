! this subroutine finds an adaptive grid cell corresponding to an input position
! this subroutine expects a position out of bounds of the dyncells
!
! INPUT         pos             position
! OUTPUT        actual_cell     returns an index of the current cell
SUBROUTINE find_dyn_cell1(pos, actual_cell)

USE types
USE constants
IMPLICIT NONE

! input variables
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: pos
! subroutine returns number of the actual cell
INTEGER                                 :: actual_cell
! basic cell variables
INTEGER, DIMENSION(const_dimofspace)                   :: bcell
INTEGER                                 :: bindex
! local actual cell
INTEGER                                 :: actCell
! parameters of subcells of dyngrid ijk
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: subcells_width
INTEGER                                         :: subind_x, subind_y, subind_z
INTEGER                                         :: sub_nx, sub_ny, sub_nz
DOUBLE PRECISION                                :: rat1, rat2, rat3

INTEGER, DIMENSION(const_dimofspace)                           :: n_cell

DOUBLE PRECISION, PARAMETER                     :: epsilon0 = 1e-6
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pobcw

n_cell = (/ nx_cell, ny_cell, nz_cell /)

! DO ind_I = 1,3
!  write(*,*) 'find_dyn_cell1: pos  = ', pos(ind_I)
!  IF(abs(pos(ind_I)) < epsilon0) THEN
!   pos(ind_I) = 0.D0
!  END IF
! END DO

! bounds test

IF(pos(ind_x) > xmax .or. pos(ind_x) < -xmax) THEN
 actual_cell = -99
 RETURN
END IF
IF(pos(ind_y) > ymax .or. pos(ind_y) < -ymax) THEN
 actual_cell = -99
 RETURN
END IF
IF(pos(ind_z) > zmax .or. pos(ind_z) < -zmax) THEN
 actual_cell = -99
 RETURN
END IF

! firstly we can compute which basic cell this point contains
pobcw = pos(:)/basic_cell_width(:)
bcell(ind_x) = FLOOR(pobcw(ind_x) + dble(nx_cell)/2.D0) + 1
bcell(ind_y) = FLOOR(pobcw(ind_y) + dble(ny_cell)/2.D0) + 1
bcell(ind_z) = FLOOR(pobcw(ind_z) + dble(nz_cell)/2.D0) + 1

! write(*,*) 'find_dyn_cell1: bcell = ', bcell, ' pobcw = ', pobcw
! write(*,*) 'find_dyn_cell1: pos(1)/width(1) = ', pos(:)/basic_cell_width(:), ' n_x/2 = ', dble(nx_cell)/2.D0
! write(*,*) 'find_dyn_cell1: (...) = ', pos(1)/basic_cell_width(1) - FLOOR(pos(1)/basic_cell_width(1))

! correction for the boundaries
IF(pos(ind_x) == xmax) THEN
 bcell(ind_x) = bcell(ind_x) - 1
END IF
IF(pos(ind_y) == ymax) THEN
 bcell(ind_y) = bcell(ind_y) - 1
END IF
IF(pos(ind_z) == zmax) THEN
 bcell(ind_z) = bcell(ind_z) - 1
END IF

IF(pos(ind_x) < xmax .and. bcell(ind_x) == nx_cell + 1) THEN
 bcell(ind_x) = bcell(ind_x) - 1
 write(*,*) 'find_dyn_cell1: korekce na bunku ve smeru x'
END IF
IF(pos(ind_y) < xmax .and. bcell(ind_y) == nx_cell + 1) THEN
 bcell(ind_y) = bcell(ind_y) - 1
 write(*,*) 'find_dyn_cell1: korekce na bunku ve smeru y'
END IF
IF(pos(ind_z) < xmax .and. bcell(ind_z) == nx_cell + 1) THEN
 bcell(ind_z) = bcell(ind_z) - 1
 write(*,*) 'find_dyn_cell1: korekce na bunku ve smeru z'
END IF

! write(*,*) 'find_dyn_cell1: bcell = ', bcell

! index of the given basic cell
bindex = (bcell(ind_x) - 1) * ny_cell * nz_cell + (bcell(ind_y) - 1) * nz_cell + bcell(ind_z)
! initial setting of the local variable corresponding to the actual cell
actCell = bindex
! if there is no dynamical cell in the given basic cell

IF(dyn_cell(actCell)%up_cell == 0) THEN
 actual_cell = actCell
 RETURN 
! we will move upper in the dyncell tree otherwise
END IF
SELECT CASE(dyngrid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dyngrid 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(1)
! we are looking for the given cell in dyncell tree
actCell = dyn_cell(actCell)%up_cell
DO
 ! did we found the given cell containing the given point?
 IF((pos(ind_x) .GE. dyn_cell(actCell)%corner(ind_x)) .AND. &
  (pos(ind_x) .LE. dyn_cell(actCell)%corner(ind_x) + dyn_cell(actCell)%width(ind_x)) .AND. &
    (pos(ind_y) .GE. dyn_cell(actCell)%corner(ind_y)) .AND. &
    (pos(ind_y) .LE. dyn_cell(actCell)%corner(ind_y) + dyn_cell(actCell)%width(ind_y)) .AND. &
    (pos(ind_z) .GE. dyn_cell(actCell)%corner(ind_z)) .AND. &
    (pos(ind_z) .LE. dyn_cell(actCell)%corner(ind_z) + dyn_cell(actCell)%width(ind_z))) THEN
  ! we have found a cell containing the given point
  ! is this cell on the top of the dyncell tree?
  IF(dyn_cell(actCell)%up_cell == 0) THEN
   actual_cell = actCell
   EXIT
  ! we have to move to the higher level of the dyncell tree
  ELSE
   actCell = dyn_cell(actCell)%up_cell
  END IF
 ! we are not in the right cell, we have to move to the next cell
 ! in the dynamical cells tree
 ELSE
  IF(actCell + 1 .LE. dyn_cell(dyn_cell(actCell)%down_cell)%up_cell + 7) THEN
   actCell = actCell + 1
  ELSE
   STOP 'error in the next dynamical cell calculating'
  END IF
 END IF
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dyngrid ijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 subcells_width = dyn_cell(dyn_cell(actCell)%up_cell)%width
 ! write(*,*) 'find_dyn_cell1: dyn_cell(actCell)%width / subcells_width = ', dyn_cell(actCell)%width(1) / subcells_width(1), &
 ! dyn_cell(actCell)%width(2) / subcells_width(2), dyn_cell(actCell)%width(3) / subcells_width(3)
! write(*,*) 'find_dyn_cell1: ROUND(dyn_cell(actCell)%width(1) / subcells_width(1)) = ', &
!  dyn_cell(actCell)%width(1) / subcells_width(1)
!  dyn_cell(actCell)%width(2) / subcells_width(2)
!  dyn_cell(actCell)%width(3) / subcells_width(3)
 rat1 = dyn_cell(actCell)%width(ind_x) / subcells_width(ind_x)
 rat2 = dyn_cell(actCell)%width(ind_y) / subcells_width(ind_y)
 rat3 = dyn_cell(actCell)%width(ind_z) / subcells_width(ind_z)
!  write(*,*) 'find_dyn_cell1: rat1 = ', rat1, ' rat2 = ', rat2, ' rat3 = ', rat3
 IF(MODULO(rat1,1.0) > 0.5) THEN
  sub_nx = CEILING(rat1)
 ELSE IF(MODULO(rat1,1.0) <= 0.5 .AND. MODULO(rat1,1.0) /= 0.0) THEN
  sub_nx = FLOOR(rat1)
 ELSE IF(MODULO(rat1,1.0) == 0.0) THEN
  sub_nx = INT(rat1)
 END IF
 IF(MODULO(rat2,1.0) > 0.5) THEN
  sub_ny = CEILING(rat2)
 ELSE IF(MODULO(rat2,1.0) <= 0.5 .AND. MODULO(rat2,1.0) /= 0.0) THEN
  sub_ny = FLOOR(rat2)
 ELSE IF(MODULO(rat2,1.0) == 0.0) THEN
  sub_ny = INT(rat2)
 END IF
 IF(MODULO(rat3,1.0) > 0.5) THEN
  sub_nz = CEILING(rat3)
 ELSE IF(MODULO(rat3,1.0) <= 0.5 .AND. MODULO(rat3,1.0) /= 0.0) THEN
  sub_nz = FLOOR(rat3)
 ELSE IF(MODULO(rat3,1.0) == 0.0) THEN
  sub_nz = INT(rat3)
 END IF

 subind_x = FLOOR((pos(ind_x) - dyn_cell(actcell)%corner(ind_x))/subcells_width(ind_x)) + 1
 subind_y = FLOOR((pos(ind_y) - dyn_cell(actcell)%corner(ind_y))/subcells_width(ind_y)) + 1
 subind_z = FLOOR((pos(ind_z) - dyn_cell(actcell)%corner(ind_z))/subcells_width(ind_z)) + 1

 actual_cell = dyn_cell(actcell)%up_cell + &
        sub_ny * sub_nz * (subind_x - 1) + &
        sub_nz * (subind_y - 1) + subind_z - 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! default case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE DEFAULT
 STOP 'find_dyn_cell1: choice of the dyngrid type is not known'
END SELECT

END SUBROUTINE find_dyn_cell1
