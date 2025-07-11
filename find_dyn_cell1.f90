! this subroutine finds an adaptive grid cell corresponding to an input position
! this subroutine expects a position out of bounds of the dyncells
!
! INPUT         pos             position
! OUTPUT        obtained_cell     returns an index of the current cell
!
! RETURN POINT: 4x
SUBROUTINE find_dyn_cell1(pos, obtained_cell)

USE types
USE constants
IMPLICIT NONE

! input variables
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: pos
! subroutine returns number of the obtained cell
INTEGER                                 :: obtained_cell
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

DOUBLE PRECISION, PARAMETER                     :: epsilon0 = 1e-15
! DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pobcw
DOUBLE PRECISION, DIMENSION(const_dimofspace)                   :: corner, width, upcorner
INTEGER                                                         :: ind_I, ind_J
DOUBLE PRECISION, PARAMETER                     :: mininum = 1e5

INTEGER                                         :: cur_cell, test_index
INTEGER                                         :: down_cell, up_cell, init_down_cell

DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_delta, cur_centre
LOGICAL                                         :: sign_x, sign_y, sign_z

n_cell = (/ nx_cell, ny_cell, nz_cell /)

! bounds test
IF(pos(ind_x) > xmax .or. pos(ind_x) < -xmax) THEN
 obtained_cell = -99
 RETURN
END IF
IF(pos(ind_y) > ymax .or. pos(ind_y) < -ymax) THEN
 obtained_cell = -99
 RETURN
END IF
IF(pos(ind_z) > zmax .or. pos(ind_z) < -zmax) THEN
 obtained_cell = -99
 RETURN
END IF

! firstly we can compute which basic cell this point contains
bcell(ind_x) = FLOOR(pos(ind_x)/basic_cell_width(ind_x) + dble(nx_cell)/2.D0 + epsilon0) + 1
bcell(ind_y) = FLOOR(pos(ind_y)/basic_cell_width(ind_y) + dble(ny_cell)/2.D0 + epsilon0) + 1
bcell(ind_z) = FLOOR(pos(ind_z)/basic_cell_width(ind_z) + dble(nz_cell)/2.D0 + epsilon0) + 1
IF(debug == 2) THEN
 write(*,*) 'find_dyn_cell1: bcell calculation'
 write(*,*) 'find_dyn_cell1: pos = ', pos
 write(*,*) 'find_dyn_cell1: pos/R_inf = ', pos/R_inf
 write(*,*) 'find_dyn_cell1: basic_cell_width = ', basic_cell_width
 write(*,*) 'find_dyn_cell1: nx_cell = ', nx_cell, ' ny_cell = ', ny_cell, ' nz_cell = ', nz_cell
 write(*,*) 'find_dyn_cell1: (pobcw(ind_x) + dble(nx_cell)/2.D0) = ', &
  (pos(ind_x)/basic_cell_width(ind_x) + dble(nx_cell)/2.D0)
 write(*,*) 'find_dyn_cell1: (pobcw(ind_y) + dble(ny_cell)/2.D0) = ', &
  (pos(ind_y)/basic_cell_width(ind_y) + dble(ny_cell)/2.D0)
 write(*,*) 'find_dyn_cell1: (pobcw(ind_z) + dble(nz_cell)/2.D0) = ', &
  (pos(ind_z)/basic_cell_width(ind_z) + dble(nz_cell)/2.D0)
 write(*,*) 'find_dyn_cell1: (pobcw(ind_x) + dble(nx_cell)/2.D0 + epsilon0) = ', &
  (pos(ind_x)/basic_cell_width(ind_x) + dble(nx_cell)/2.D0 + epsilon0)
 write(*,*) 'find_dyn_cell1: (pobcw(ind_y) + dble(ny_cell)/2.D0 + epsilon0) = ', &
  (pos(ind_y)/basic_cell_width(ind_y) + dble(ny_cell)/2.D0 + epsilon0)
 write(*,*) 'find_dyn_cell1: (pobcw(ind_z) + dble(nz_cell)/2.D0 + epsilon0) = ', &
  (pos(ind_z)/basic_cell_width(ind_z) + dble(nz_cell)/2.D0 + epsilon0)
 write(*,*) 'find_dyn_cell1: bcell = ', bcell
 write(*,*) 'find_dyn_cell1: ___________________________________________'
END IF

IF(pos(ind_x) <= xmax .and. bcell(ind_x) == nx_cell + 1) THEN
 bcell(ind_x) = bcell(ind_x) - 1
END IF
IF(pos(ind_y) <= xmax .and. bcell(ind_y) == nx_cell + 1) THEN
 bcell(ind_y) = bcell(ind_y) - 1
END IF
IF(pos(ind_z) <= xmax .and. bcell(ind_z) == nx_cell + 1) THEN
 bcell(ind_z) = bcell(ind_z) - 1
END IF

! index of the given basic cell
bindex = (bcell(ind_x) - 1) * ny_cell * nz_cell + (bcell(ind_y) - 1) * nz_cell + bcell(ind_z)
! initial setting of the local variable corresponding to the actual cell
actCell = bindex
! if there is no dynamical cell in the given basic cell
IF(actCell > nx_cell * ny_cell * nz_cell) THEN
 write(*,*) 'find_dyn_cell1: xmax = ', xmax, ' ymax = ', ymax, ' zmax = ', zmax
 write(*,*) 'find_dyn_cell1: w_x = ', basic_cell_width(ind_x), ' w_y  = ', basic_cell_width(ind_y), &
  ' w_z = ', basic_cell_width(ind_z)
 write(*,*) 'find_dyn_cell1: x/w_x', pos(ind_x)/basic_cell_width(ind_x)
 write(*,*) 'find_dyn_cell1: y/w_y', pos(ind_y)/basic_cell_width(ind_y)
 write(*,*) 'find_dyn_cell1: z/w_z', pos(ind_z)/basic_cell_width(ind_z)
 write(*,*) 'find_dyn_cell1: pos_x/xmax = ', pos(ind_x)/xmax, ' pos_y/ymax = ', pos(ind_y)/ymax, &
  ' pos_z/zmax = ', pos(ind_z)/zmax
 write(*,*) 'find_dyn_cell1: bcell_x = ', bcell(ind_x), ' bcell_y = ', bcell(ind_y), &
  ' bcell_z = ', bcell(ind_z)
END IF

IF(debug == 2) THEN
 corner = dyn_cell(actCell)%corner
 upcorner = dyn_cell(actCell)%upcorner
 width = dyn_cell(actCell)%width
 DO ind_I = 1, const_dimofspace
  IF(((pos(ind_I) < corner(ind_I) - mininum) .OR. (pos(ind_I) > upcorner(ind_I)+ mininum))) THEN
   write(*,*) 'find_dyn_cell1: test #1'
   write(*,*) 'find_dyn_cell1: cell starting = ', corner
   write(*,*) 'find_dyn_cell1: packet pos = ', pos
   write(*,*) 'find_dyn_cell1: cell ending = ', upcorner
   write(*,*) 'find_dyn_cell1: ind_I = ', ind_I
   write(*,*) 'find_dyn_cell1: packet is not located inside the propagation cell'
   STOP
  END IF
 END DO
END IF

IF(dyn_cell(actCell)%up_cell == 0) THEN
 obtained_cell = actCell
 IF(debug == 2) THEN
  corner = dyn_cell(obtained_cell)%corner
  upcorner = dyn_cell(obtained_cell)%upcorner
  DO ind_I = 1, const_dimofspace
   IF(((pos(ind_I) < corner(ind_I) - mininum) .OR. (pos(ind_I) > upcorner(ind_I) + mininum))) THEN
    write(*,*) 'find_dyn_cell1: test #2'
    write(*,*) 'find_dyn_cell1: cell starting = ', corner
    write(*,*) 'find_dyn_cell1: packet pos = ', pos
    write(*,*) 'find_dyn_cell1: cell ending = ', upcorner
    write(*,*) 'find_dyn_cell1: ind_I = ', ind_I
    write(*,*)  'find_dyn_cell1: packet is not located inside the propagation cell'
    STOP
   END IF
  END DO
 END IF
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
! we are looking for the given cell in the dyncell tree
ind_J = 0
cur_centre = (dyn_cell(actCell)%upcorner + dyn_cell(actCell)%corner) / 2.D0

cur_cell = actCell
init_down_cell = actCell

DO
 ! we have to move to the higher level of the dyncell tree
 ind_J = ind_J + 1
 write(*,*) 'find_dyn_cell1: loop ind_J = ', ind_J, ' cur_cell = ', cur_cell
 ! IF(cur_cell == dyn_cell(cur_cell)%up_cell) THEN
 !  write(*,*) 'find_dyn_cell1: ERROR upper cell == actCell'
 !  STOP 'find_dyn_cell1'
 ! END IF
 ! did we found the given cell containing the given point?
 ! DO ind_I = 0,7
 !  test_index = dyn_cell(cur_cell)%up_cell + ind_I
 !  corner = dyn_cell(test_index)%corner
 !  upcorner = dyn_cell(test_index)%upcorner
 !  write(*,*) 'find_dyn_cell1: ind_I + 1 = ', ind_I + 1
 !  write(*,*) 'find_dyn_cell1: corner = ', corner
 !  write(*,*) 'find_dyn_cell1: pos = ', pos
 !  write(*,*) 'find_dyn_cell1: upcorner = ', upcorner
 ! END DO

 IF((pos(ind_x) .GE. dyn_cell(cur_cell)%corner(ind_x) - mininum) .AND. &
    (pos(ind_x) .LE. dyn_cell(cur_cell)%upcorner(ind_x) + mininum) .AND. &
    (pos(ind_y) .GE. dyn_cell(cur_cell)%corner(ind_y) - mininum) .AND. &
    (pos(ind_y) .LE. dyn_cell(cur_cell)%upcorner(ind_y) + mininum) .AND. &
    (pos(ind_z) .GE. dyn_cell(cur_cell)%corner(ind_z) - mininum) .AND. &
    (pos(ind_z) .LE. dyn_cell(cur_cell)%upcorner(ind_z) + mininum)) THEN
  ! we have found a cell containing the given point
  ! is this cell on the top of the dyncell tree?
  IF(dyn_cell(cur_cell)%up_cell == 0) THEN
   obtained_cell = cur_cell
   write(*,*) 'find_dyn_cell1: obtained_cell = ', obtained_cell
   EXIT
  ! we have to move to the higher level of the dyncell tree
  ELSE
   init_down_cell = dyn_cell(cur_cell)%up_cell
   cur_cell = dyn_cell(cur_cell)%up_cell
  END IF
 ! we are not in the right cell, we have to move to the next cell
 ! in the dynamical cells tree
 ELSE
  IF(cur_cell + 1 .LE. init_down_cell + 7) THEN
   cur_cell = cur_cell + 1
  ELSE
   ! write(*,*) 'find_dyn_cell1: up_cell = ', dyn_cell(init_down_cell)%up_cell
   ! down_cell = dyn_cell(cur_cell)%down_cell
   ! write(*,*) 'find_dyn_cell1: low cell corner = ', dyn_cell(down_cell)%corner
   ! up_cell = dyn_cell(down_cell)%up_cell
   ! DO ind_I = 1,8
   !  write(*,*) 'find_dyn_cell1: corner for pgi_index = ', up_cell, &
   !   ' corner = ', (dyn_cell(up_cell)%corner)
   ! up_cell = up_cell + 1
   ! END DO
   DO ind_I = 0,7
    test_index = dyn_cell(init_down_cell)%up_cell + ind_I
    corner = dyn_cell(test_index)%corner
    upcorner = dyn_cell(test_index)%upcorner
    write(*,*) 'find_dyn_cell1: ind_I + 1 = ', ind_I + 1
    write(*,*) 'find_dyn_cell1: corner = ', corner
    write(*,*) 'find_dyn_cell1: pos = ', pos
    write(*,*) 'find_dyn_cell1: upcorner = ', upcorner
   END DO
   STOP 'error in the next dynamical cell calculating'
  END IF
 END IF

 ! we have found a cell containing the given point
 ! is this cell on the top of the dyncell tree?
 ! we are not in the right cell, we have to move to the next cell
 ! in the dynamical cells tree
 write(*,*) 'find_dyn_cell1: cur_cell = ', cur_cell
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

 obtained_cell = dyn_cell(actcell)%up_cell + &
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
