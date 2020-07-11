SUBROUTINE next_cell_up(pack_index, dist, cell_down, next_cell)

USE types

IMPLICIT NONE

! input variables
DOUBLE PRECISION                        :: dist
INTEGER                                 :: cell_down, pack_index
! outpu variables
INTEGER                                 :: next_cell
! local variable next cell
! cross the surface
DOUBLE PRECISION, DIMENSION(3)          :: cross_pos
INTEGER                                 :: cross
! actual cell
INTEGER                                 :: act_cell
DOUBLE PRECISION, DIMENSION(3)          :: corner, width
INTEGER                                 :: upper_cell
! parameters of subcells of dyngrid ijk
DOUBLE PRECISION, DIMENSION(3)          :: subcells_width
INTEGER                                 :: subind_x, subind_y, subind_z
INTEGER                                         :: sub_nx, sub_ny, sub_nz
DOUBLE PRECISION                                :: rat1, rat2, rat3

! IF(package(pack_index)%pos > R_inf .OR. package(pack_index)%pos < R_star) THEN
!  next_cell = 

act_cell = cell_down
cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist
cross = package(pack_index)%next_cross
! is the photon close to the edge of the propagation grid?
IF(act_cell < 0) THEN
 next_cell = act_cell
 RETURN
END IF
IF(act_cell > SIZE(dyn_cell)) THEN
 write(*,*) 'next_cell_up: act_cell = ', act_cell
 CALL abort()
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
 width = dyn_cell(act_cell)%width
 upper_cell = dyn_cell(act_cell)%up_cell
! print*, 'next_cell_up: act_cell = ', act_cell, ' upper_cell = ', upper_cell, 'cell_down = ', cell_down
 IF(upper_cell == 0) THEN
  next_cell = act_cell
  EXIT
 ELSE IF(upper_cell > 0) THEN
  ! we have to find which cell in the higher level corresponds to the cross point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in x direction
   IF(cross == posx .OR. cross == negx) THEN
 !   print*, 'next_cell_up: cross = ', cross
    ! lower front cell
    IF(cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0 .AND. &
       cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
 !   print*, 'next_cell_up: cross 1'
     if(cross == posx) act_cell = upper_cell
     if(cross == negx) act_cell = upper_cell + 1
    ! lower rear cell
    ELSE IF(cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2) .AND. &
       cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
 !   print*, 'next_cell_up: cross 2'
     if(cross == posx) act_cell = upper_cell + 2
     if(cross == negx) act_cell = upper_cell + 3
    ! upper front cell
    ELSE IF(cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0 .AND. &
       cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
 !   print*, 'next_cell_up: cross 3'
     if(cross == posx) act_cell = upper_cell + 4
     if(cross == negx) act_cell = upper_cell + 5
    ! upper rear cell
    ELSE IF(cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2) .AND. &
       cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
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
   IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
!   print*, 'next_cell_up: cross 1'
    if(cross == posy) act_cell = upper_cell
    if(cross == negy) act_cell = upper_cell + 2
   ! lower right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= corner(3) + width(3) / 2.D0) THEN
!   print*, 'next_cell_up: cross 2'
    if(cross == posy) act_cell = upper_cell + 1
    if(cross == negy) act_cell = upper_cell + 3
   ! upper left cell
   ELSE IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
!   print*, 'next_cell_up: cross 3'
    if(cross == posy) act_cell = upper_cell + 4
    if(cross == negy) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(3) >= corner(3) + width(3) / 2.D0 .AND. cross_pos(3) <= corner(3) + width(3)) THEN
!   print*, 'next_cell_up: cross 4'
    if(cross == posy) act_cell = upper_cell + 5
    if(cross == negy) act_cell = upper_cell + 7
   ELSE
    write(*,*) 'next_cell_up: no cell was found'
    write(*,*) 'pos = ', norm2(package(pack_index)%pos) / R_inf
    CALL abort()
   END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in z direction
  ELSE IF(cross == posz .OR. cross == negz) THEN
!   print*, 'next_cell_up: cross = ', cross
   ! front left cell
   IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0) THEN
!   print*, 'next_cell_up: cross 1'
    if(cross == posz) act_cell = upper_cell
    if(cross == negz) act_cell = upper_cell + 4
   ! lower right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(2) >= corner(2) .AND. cross_pos(2) <= corner(2) + width(2) / 2.D0) THEN
!   print*, 'next_cell_up: cross 2'
    if(cross == posz) act_cell = upper_cell + 1
    if(cross == negz) act_cell = upper_cell + 5
   ! upper left cell
   ELSE IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= corner(1) + width(1) / 2.D0 .AND. &
      cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2)) THEN
!   print*, 'next_cell_up: cross 3'
    if(cross == posz) act_cell = upper_cell + 2
    if(cross == negz) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(1) >= corner(1) + width(1) / 2.D0 .AND. cross_pos(1) <= corner(1) + width(1) .AND. &
      cross_pos(2) >= corner(2) + width(2) / 2.D0 .AND. cross_pos(2) <= corner(2) + width(2)) THEN
!   print*, 'next_cell_up: cross 4'
    if(cross == posz) act_cell = upper_cell + 3
    if(cross == negz) act_cell = upper_cell + 7
   ELSE
    write(*,*) 'next_cell_up: no cell was found'
    write(*,*) 'pos = ', norm2(package(pack_index)%pos) / R_star
    write(*,*) 'cor_z = ', corner(3)/R_star, ' pos_z = ', package(pack_index)%pos(3)/R_star, ' cor_z+w = ', (corner(3) + &
    width(3))/R_star
    CALL abort()
   END IF
  END IF
 END IF  
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
 subcells_width = dyn_cell(dyn_cell(act_cell)%up_cell)%width
 subcells_width = dyn_cell(dyn_cell(act_cell)%up_cell)%width
 rat1 = dyn_cell(act_cell)%width(1) / subcells_width(1)
 rat2 = dyn_cell(act_cell)%width(2) / subcells_width(2)
 rat3 = dyn_cell(act_cell)%width(3) / subcells_width(3)
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

 subind_x = FLOOR((cross_pos(1) - dyn_cell(act_cell)%corner(1))/subcells_width(1)) + 1
 subind_y = FLOOR((cross_pos(2) - dyn_cell(act_cell)%corner(2))/subcells_width(2)) + 1
 subind_z = FLOOR((cross_pos(3) - dyn_cell(act_cell)%corner(3))/subcells_width(3)) + 1
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
! print*, 'next_cell_up:', subind_x, subind_y, subind_z
 next_cell = dyn_cell(act_cell)%up_cell + &
        sub_ny * sub_nz * (subind_x - 1) + &
        sub_nz * (subind_y - 1) + subind_z - 1
! this case occurs also when the SBR resonance distance is
! looking for another boundary
IF(next_cell > SIZE(dyn_cell)) THEN
 next_cell = -99
END IF
! print*, 'next_cell_up: next_cell = ', next_cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! default case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE DEFAULT
 STOP 'next_cell_up: wrong choice of a dynamical grid type'
END SELECT

END SUBROUTINE next_cell_up
