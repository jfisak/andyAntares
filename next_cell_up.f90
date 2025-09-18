! finds the next cell in the upper levels of the adaptive propGrid
!
! INPUT: pack_index(INT): index of a packet
!        dist(DBLE): distance to the next propGrid cross
!        cell_down(INT): index of the next propGrid cell in the lower levels
! OUTPUT: next_cell(INT): calculated next propGrid cell index
!
! RETURN point 2X
!
SUBROUTINE next_cell_up(pack_index, dist, cell_down, next_cell)

USE types
USE constants

IMPLICIT NONE

! input variables
DOUBLE PRECISION                        :: dist
INTEGER                                 :: cell_down, pack_index
! outpu variables
INTEGER                                 :: next_cell
! local variable next cell
! cross the surface
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cross_pos
INTEGER                                 :: cross
! actual cell
INTEGER                                 :: act_cell, ind_I
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: corner, upcorner, centre
INTEGER                                 :: upper_cell
! parameters of subcells of dyngrid ijk
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: subcells_width
INTEGER                                 :: subind_x, subind_y, subind_z
INTEGER                                         :: sub_nx, sub_ny, sub_nz

INTEGER, DIMENSION(const_dimofspace)                   :: bcell
INTEGER                                 :: bindex
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: pos

DOUBLE PRECISION, PARAMETER                     :: epsilon0 = 1e-6
! IF(package(pack_index)%pos > R_inf .OR. package(pack_index)%pos < R_star) THEN
!  next_cell = 

act_cell = cell_down
cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist
cross = package(pack_index)%next_cross
IF(debug == 2) THEN
 write(*,*) 'next_cell_up: __________begin______________________________'
 write(*,*) 'next_cell_up: cross = ', cross
 write(*,*) 'next_cell_up: pos = ', package(pack_index)%pos
 write(*,*) 'next_cell_up: neighbors = ', dyn_cell(act_cell)%neighbor
 write(*,*) 'next_cell_up: up cell = ', dyn_cell(act_cell)%up_cell
 write(*,*) 'next_cell_up: __________begin______________________________'
END IF
! is the photon close to the edge of the propagation grid?
IF(act_cell < 0) THEN
 next_cell = act_cell
 RETURN
END IF

corner = dyn_cell(act_cell)%corner
upcorner = dyn_cell(act_cell)%upcorner
IF(cross == posx) THEN
 IF(cross_pos(ind_x) /= corner(ind_x)) THEN
  IF(debug == 2) write(*,*) 'next_cell_up: correcting position for cross = ', cross, ' ind ', ind_x
  cross_pos(ind_x) = corner(ind_x)
 END IF
ELSE IF(cross == negx) THEN
 IF(cross_pos(ind_x) /= upcorner(ind_x)) THEN
  IF(debug == 2) write(*,*) 'next_cell_up: correcting position for cross = ', cross, ' ind ', ind_x
  cross_pos(ind_x) = upcorner(ind_x)
 END IF
ELSE IF(cross == posy) THEN
 IF(cross_pos(ind_y) /= corner(ind_y)) THEN
  IF(debug == 2) write(*,*) 'next_cell_up: correcting position for cross = ', cross, ' ind ', ind_y
  cross_pos(ind_y) = corner(ind_y)
 END IF
ELSE IF(cross == negy) THEN
 IF(cross_pos(ind_y) /= upcorner(ind_y)) THEN
  IF(debug == 2) write(*,*) 'next_cell_up: correcting position for cross = ', cross, ' ind ', ind_y
  cross_pos(ind_y) = upcorner(ind_y)
 END IF
ELSE IF(cross == posz) THEN
 IF(cross_pos(ind_z) /= corner(ind_z)) THEN
  IF(debug == 2) write(*,*) 'next_cell_up: correcting position for cross = ', cross, ' ind ', ind_z
  cross_pos(ind_z) = corner(ind_z)
 END IF
ELSE IF(cross == negz) THEN
 IF(cross_pos(ind_z) /= upcorner(ind_z)) THEN
  IF(debug == 2) write(*,*) 'next_cell_up: correcting position for cross = ', cross, ' ind ', ind_z
  cross_pos(ind_z) = upcorner(ind_z)
 END IF
END IF
IF(act_cell > SIZE(dyn_cell)) THEN
 ! write(*,*) 'next_cell_up: act_cell = ', act_cell
 ! CALL abort()
END IF

SELECT CASE(dyngrid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dynamical grid type 8
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! there is a bug in this type of cells
CASE(1)
!print*, 'next_cell_up: act_cell = ', act_cell, ' cross_pos = ', cross_pos, 'cross = ', cross
DO
 corner = dyn_cell(act_cell)%corner
 upcorner = dyn_cell(act_cell)%upcorner
 upper_cell = dyn_cell(act_cell)%up_cell
 IF(upper_cell > 0) centre = dyn_cell(upper_cell)%upcorner
 ! write(*,*) 'next_cell_up: corner = ', corner
 IF(upper_cell > 0) write(*,*) 'next_cell_up: centre = ', centre
 ! write(*,*) 'next_cell_up: cross_pos = ', cross_pos
 ! write(*,*) 'next_cell_up: upcorner = ', upcorner
 ! write(*,*) 'next_cell_up: upper_cell = ', upper_cell
 IF(upper_cell == 0) THEN
  next_cell = act_cell
  IF(debug == 2) THEN
   corner = dyn_cell(next_cell)%corner
   upcorner = dyn_cell(next_cell)%upcorner
   ! write(*,*) '_____calculated cell__________________'
   ! write(*,*) 'next_cell_up: next_cell = ', next_cell
   ! write(*,*) 'next_cell_up: corner = ', corner
   ! write(*,*) 'next_cell_up: cross_pos = ', cross_pos
   ! write(*,*) 'next_cell_up: upcorner = ', upcorner
   ! write(*,*) '______________________________________'
   DO ind_I = 1, const_dimofspace
    IF(((cross_pos(ind_I) < corner(ind_I) ) .OR. (cross_pos(ind_I) > upcorner(ind_I) ))) THEN
     upper_cell = dyn_cell(act_cell)%up_cell
     IF(upper_cell > 0) centre = dyn_cell(upper_cell)%upcorner
     write(*,*) 'next_cell_up: _________________________________________________'
     IF(upper_cell > 0) write(*,*) 'next_cell_up: centre = ', centre
     write(*,*) 'next_cell_up: cell starting = ', corner
     write(*,*) 'next_cell_up: packet pos = ', cross_pos
     write(*,*) 'next_cell_up: cell ending = ', upcorner
     write(*,*) 'next_cell_up: _________________________________________________'
     write(*,*) 'I = ', ind_I
     STOP 'next_cell_up: packet is not located inside the propagation cell'
    END IF
   END DO
  END IF
  EXIT
 ELSE IF(upper_cell > 0) THEN
  write(*,*) 'next_cell_up: cross = ', cross
  ! we have to find which cell in the higher level corresponds to the cross point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in x direction
   IF(cross == posx .OR. cross == negx) THEN
 !   print*, 'next_cell_up: cross = ', cross
    ! lower front cell
    IF(cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= centre(ind_y)  .AND. &
       cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= centre(ind_z) ) THEN
 !   print*, 'next_cell_up: cross 1'
     if(cross == posx) act_cell = upper_cell
     if(cross == negx) act_cell = upper_cell + 1
    ! lower rear cell
    ELSE IF(cross_pos(ind_y) >= centre(ind_y) .AND. cross_pos(ind_y) <= upcorner(ind_y) .AND. &
            cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= centre(ind_z)) THEN
 !   print*, 'next_cell_up: cross 2'
     if(cross == posx) act_cell = upper_cell + 2
     if(cross == negx) act_cell = upper_cell + 3
    ! upper front cell
    ELSE IF(cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= centre(ind_y) .AND. &
       cross_pos(ind_z) >= centre(ind_z) .AND. cross_pos(ind_z) <= upcorner(ind_z)) THEN
 !   print*, 'next_cell_up: cross 3'
     if(cross == posx) act_cell = upper_cell + 4
     if(cross == negx) act_cell = upper_cell + 5
    ! upper rear cell
    ELSE IF(cross_pos(ind_y) >= centre(ind_y) .AND. cross_pos(ind_y) <= upcorner(ind_y) .AND. &
       cross_pos(ind_z) >= centre(ind_z) .AND. cross_pos(ind_z) <= upcorner(ind_z)) THEN
 !   print*, 'next_cell_up: cross 4'
     if(cross == posx) act_cell = upper_cell + 6
     if(cross == negx) act_cell = upper_cell + 7
    ELSE
     write(*,*) 'next_cell_up: no cell was found'
     CALL abort()
    END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in y direction
  ELSE IF(cross == posy .OR. cross == negy) THEN
!   print*, 'next_cell_up: cross = ', cross
   ! lower left cell
   IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= centre(ind_x)  .AND. &
      cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= centre(ind_z) ) THEN
!   print*, 'next_cell_up: cross 1'
    if(cross == posy) act_cell = upper_cell
    if(cross == negy) act_cell = upper_cell + 2
   ! lower right cell
   ELSE IF(cross_pos(ind_x) >= centre(ind_x) .AND. cross_pos(ind_x) <= upcorner(ind_x) .AND. &
      cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= centre(ind_z)) THEN
!   print*, 'next_cell_up: cross 2'
    if(cross == posy) act_cell = upper_cell + 1
    if(cross == negy) act_cell = upper_cell + 3
   ! upper left cell
   ELSE IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= centre(ind_x) .AND. &
      cross_pos(ind_z) >= centre(ind_z) .AND. cross_pos(ind_z) <= upcorner(ind_z)) THEN
!   print*, 'next_cell_up: cross 3'
    if(cross == posy) act_cell = upper_cell + 4
    if(cross == negy) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(ind_x) >= centre(ind_x) .AND. cross_pos(ind_x) <= upcorner(ind_x) .AND. &
      cross_pos(ind_z) >= centre(ind_z) .AND. cross_pos(ind_z) <= upcorner(ind_z)) THEN
!   print*, 'next_cell_up: cross 4'
    if(cross == posy) act_cell = upper_cell + 5
    if(cross == negy) act_cell = upper_cell + 7
   ELSE
    write(*,*) 'next_cell_up: no cell was found'
    CALL abort()
   END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in z direction
  ELSE IF(cross == posz .OR. cross == negz) THEN
!   print*, 'next_cell_up: cross = ', cross
   ! front left cell
   IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= centre(ind_x) .AND. &
      cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= centre(ind_y)) THEN
!   print*, 'next_cell_up: cross 1'
    if(cross == posz) act_cell = upper_cell
    if(cross == negz) act_cell = upper_cell + 4
   ! lower right cell
   ELSE IF(cross_pos(ind_x) >= centre(ind_x) .AND. cross_pos(ind_x) <= upcorner(ind_x) .AND. &
      cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= centre(ind_y)) THEN
!   print*, 'next_cell_up: cross 2'
    if(cross == posz) act_cell = upper_cell + 1
    if(cross == negz) act_cell = upper_cell + 5
   ! upper left cell
   ELSE IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= centre(ind_x) .AND. &
      cross_pos(ind_y) >= centre(ind_y) .AND. cross_pos(ind_y) <= upcorner(ind_y)) THEN
!   print*, 'next_cell_up: cross 3'
    if(cross == posz) act_cell = upper_cell + 2
    if(cross == negz) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(ind_x) >= centre(ind_x) .AND. cross_pos(ind_x) <= upcorner(ind_x) .AND. &
      cross_pos(ind_y) >= centre(ind_y) .AND. cross_pos(ind_y) <= upcorner(ind_y)) THEN
!   print*, 'next_cell_up: cross 4'
    if(cross == posz) act_cell = upper_cell + 3
    if(cross == negz) act_cell = upper_cell + 7
   ELSE
    write(*,*) 'next_cell_up: no cell was found'
    CALL abort()
   END IF
  END IF
 END IF  
 ! write(*,*) 'next_cell_up: _________________________________________________'
 ! write(*,*) 'next_cell_up: act_cell = ', act_cell
 ! write(*,*) 'next_cell_up: _________________________________________________'
 corner = dyn_cell(act_cell)%corner
 upcorner = dyn_cell(act_cell)%upcorner
 DO ind_I = 1, const_dimofspace
  IF(((cross_pos(ind_I) < corner(ind_I) ) .OR. (cross_pos(ind_I) > upcorner(ind_I) ))) THEN
   upper_cell = dyn_cell(act_cell)%up_cell
   IF(upper_cell > 0) centre = dyn_cell(upper_cell)%upcorner
   write(*,*) 'next_cell_up: _________________________________________________'
   IF(upper_cell > 0) write(*,*) 'next_cell_up: centre = ', centre
   write(*,*) 'next_cell_up: cell starting = ', corner
   write(*,*) 'next_cell_up: packet pos = ', cross_pos
   write(*,*) 'next_cell_up: cell ending = ', upcorner
   write(*,*) 'next_cell_up: _________________________________________________'
   write(*,*) 'I = ', ind_I
   STOP 'next_cell_up: packet is not located inside the propagation cell'
  END IF
 END DO

END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! dynamical grid type ijk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 upper_cell = dyn_cell(act_cell)%up_cell
 IF(upper_cell == 0) THEN
  next_cell = act_cell
  RETURN
 END IF
 subcells_width = dyn_cell(upper_cell)%width
 subcells_width = dyn_cell(upper_cell)%width
 sub_nx = dyn_cell(act_cell)%n_sbgr(ind_x)
 sub_ny = dyn_cell(act_cell)%n_sbgr(ind_y)
 sub_nz = dyn_cell(act_cell)%n_sbgr(ind_z)

 subind_x = FLOOR((cross_pos(ind_x) - dyn_cell(act_cell)%corner(ind_x))/subcells_width(ind_x)) + 1
 subind_y = FLOOR((cross_pos(ind_y) - dyn_cell(act_cell)%corner(ind_y))/subcells_width(ind_y)) + 1
 subind_z = FLOOR((cross_pos(ind_z) - dyn_cell(act_cell)%corner(ind_z))/subcells_width(ind_z)) + 1

 pos = package(pack_index)%pos
 bcell(ind_x) = FLOOR(pos(ind_x)/basic_cell_width(ind_x) + dble(nx_cell)/2.D0) + 1
 bcell(ind_y) = FLOOR(pos(ind_y)/basic_cell_width(ind_y) + dble(ny_cell)/2.D0) + 1
 bcell(ind_z) = FLOOR(pos(ind_z)/basic_cell_width(ind_z) + dble(nz_cell)/2.D0) + 1
  
 ! index of the given basic cell
 bindex = (bcell(ind_x) - 1) * ny_cell * nz_cell + (bcell(ind_y) - 1) * nz_cell + bcell(ind_z)
 ! write(*,*) 'lower_grid = ', dyn_cell(upper_cell)%down_cell
 ! write(*,*) 'act_cell = ', act_cell, ' bindex = ', bindex

!  IF(bindex /= act_cell) THEN
!   write(*,*) 'next_cell_up: bindex = ', bindex, ' act_cell = ', act_cell
!   write(*,*) 'next_cell_up: pos = ', pos
!   write(*,*) 'next_cell_up: corner BC= ', dyn_cell(bindex)%corner
!   write(*,*) 'next_cell_up: corner2 BC= ', dyn_cell(bindex)%corner+dyn_cell(bindex)%width
!   write(*,*) 'next_cell_up: corner: ', dyn_cell(act_cell)%corner
!   write(*,*) 'next_cell_up: corner2: ', dyn_cell(act_cell)%corner+dyn_cell(act_cell)%width
!   STOP 'basic cell index do not agree with the current cell index'
!  END IF

!  IF(subind_x < 1 .or. subind_y < 1 .or. subind_z < 1 &
!   .or. subind_x > sub_nx .or. subind_y > sub_ny .or. subind_z > sub_nz) THEN
!   ! firstly we can compute which basic cell this point contains
!   write(*,*) 'next_cell_up: sub_nx = ', sub_nx, ' sub_ny = ', sub_ny, ' sub_nz = ', sub_nz
!   write(*,*) 'next_cell_up: subind_x = ', subind_x, ' subind_y = ', subind_y, ' subind_z = ', subind_z
!   write(*,*) 'lower_grid = ', dyn_cell(upper_cell)%down_cell
!   write(*,*) 'act_cell = ', act_cell, ' bindex = ', bindex
!   write(*,*) 'r/w = ', dyn_cell(act_cell)%corner(:)/subcells_width(:)
!   STOP 'subind < 1'
!  END IF

 ! write(*,*) 'next_cell_up: subind_x = ', subind_x, ' subind_y = ', subind_y, ' subind_z = ', subind_z

 IF(cross == posx) THEN
  subind_x = 1
 ELSE IF(cross == negx) THEN
  subind_x = sub_nx
 ELSE IF(cross == posy) THEN
  subind_y = 1
 ELSE IF(cross == negy) THEN
  subind_y = sub_ny
 ELSE IF(cross == posz) THEN
  subind_z = 1
 ELSE IF(cross == negz) THEN
  subind_z = sub_nz
 END IF
 next_cell = upper_cell + &
        sub_ny * sub_nz * (subind_x - 1) + &
        sub_nz * (subind_y - 1) + subind_z - 1
! this case occurs also when the SBR resonance distance is
! looking for another boundary
IF(next_cell > SIZE(dyn_cell)) THEN
 next_cell = -99
END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! default case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE DEFAULT
 STOP 'next_cell_up: wrong choice of a dynamical grid type'
END SELECT

END SUBROUTINE next_cell_up
