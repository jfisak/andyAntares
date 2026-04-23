! choose the next cell for the d-package
!
! INPUT: cur_pgi(INT): curret propGrid index
!        next_leak(INT): next direction randomly chosen
! OUTPUT: next_cell(INT): next propGrid cell
!         cross_pos(DBLE(const_dimofspace)): position of the poit in the boundary
!
! 2x RETURN POINT
!
SUBROUTINE d_choosenextcell(cur_pgi, next_leak, next_cell, cross_pos)

USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: cur_pgi, next_cell, next_leak

INTEGER                                         :: n_cell

INTEGER                                         :: ind_I

DOUBLE PRECISION                                :: rand, ran2
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: corner, width, cross_pos, upcorner

INTEGER                                         :: subind_x, subind_y, subind_z
INTEGER                                         :: sub_nx, sub_ny, sub_nz
INTEGER                                         :: act_cell, upper_cell

DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: subcells_width



corner = dyn_cell(cur_pgi)%corner
upcorner = dyn_cell(cur_pgi)%upcorner
rand = ran2(idum)
width = (corner + upcorner)/2.D0

DO ind_I = 1,3
 cross_pos(ind_I) = corner(ind_I) + rand * width(ind_I)
END DO

! write(*,*) 'd_choosenextcell: next_leak = ', next_leak

IF(next_leak == posx) THEN
 cross_pos(ind_x) = corner(ind_x) + width(ind_x)
ELSE IF(next_leak == negx) THEN
 cross_pos(ind_x) = corner(ind_x)
ELSE IF(next_leak == posy) THEN
 cross_pos(ind_y) = corner(ind_y) + width(ind_y)
ELSE IF(next_leak == negy) THEN
 cross_pos(ind_y) = corner(ind_y)
ELSE IF(next_leak == posz) THEN
 cross_pos(ind_z) = corner(ind_z) + width(ind_z)
ELSE IF(next_leak == negz) THEN
 cross_pos(ind_z) = corner(ind_z)
END IF


! next cell down
DO
 n_cell = dyn_cell(cur_pgi)%neighbor(next_leak)
 act_cell = cur_pgi
 IF(n_cell == 0) THEN
  IF(dyn_cell(act_cell)%down_cell == 0) STOP 'next_cell: no cell was found'
  act_cell = dyn_cell(act_cell)%down_cell
 ELSE 
  next_cell = n_cell
  IF(next_cell > SIZE(dyn_cell)) THEN
   CALL abort()
  END IF
  EXIT
 END IF
END DO

! if we are at the propGrid edge, we do not have to find upper cells
! RETURN POINT
IF(next_cell < 0) RETURN

! next cell up
IF(dyn_cell(next_cell)%up_cell > 0) THEN

 ! write(*,*) 'd_choosenextcell: next_cell = ', next_cell

 SELECT CASE(dyngrid)
 CASE(0)
  return
 CASE(1)
 !print*, 'next_cell_up: act_cell = ', act_cell, ' cross_pos = ', cross_pos, 'cross = ', cross
 act_cell = next_cell
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
    IF(next_leak == posx .OR. next_leak == negx) THEN
  !   print*, 'next_cell_up: next_leak = ', next_leak
     ! lower front cell
     IF(cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y) / 2.D0 .AND. &
        cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z) / 2.D0) THEN
  !   print*, 'next_cell_up: cross 1'
      if(next_leak == posx) act_cell = upper_cell
      if(next_leak == negx) act_cell = upper_cell + 1
     ! lower rear cell
     ELSE IF(cross_pos(ind_y) >= corner(ind_y) + width(ind_y) / 2.D0 .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y) .AND. &
        cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z) / 2.D0) THEN
  !   print*, 'next_cell_up: cross 2'
      if(next_leak == posx) act_cell = upper_cell + 2
      if(next_leak == negx) act_cell = upper_cell + 3
     ! upper front cell
     ELSE IF(cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y) / 2.D0 .AND. &
        cross_pos(ind_z) >= corner(ind_z) + width(ind_z) / 2.D0 .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z)) THEN
  !   print*, 'next_cell_up: cross 3'
      if(next_leak == posx) act_cell = upper_cell + 4
      if(next_leak == negx) act_cell = upper_cell + 5
     ! upper rear cell
     ELSE IF(cross_pos(ind_y) >= corner(ind_y) + width(ind_y) / 2.D0 .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y) .AND. &
        cross_pos(ind_z) >= corner(ind_z) + width(ind_z) / 2.D0 .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z)) THEN
  !   print*, 'next_cell_up: cross 4'
      if(next_leak == posx) act_cell = upper_cell + 6
      if(next_leak == negx) act_cell = upper_cell + 7
     ELSE
      ! write(*,*) 'next_cell_up: no cell was found'
      ! CALL abort()
     END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! cross in y direction
   ELSE IF(next_leak == posy .OR. next_leak == negy) THEN
 !   print*, 'next_cell_up: cross = ', cross
    ! lower left cell
    IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) / 2.D0 .AND. &
       cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z) / 2.D0) THEN
 !   print*, 'next_cell_up: cross 1'
     if(next_leak == posy) act_cell = upper_cell
     if(next_leak == negy) act_cell = upper_cell + 2
    ! lower right cell
    ELSE IF(cross_pos(ind_x) >= corner(ind_x) + width(ind_x) / 2.D0 .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) .AND. &
       cross_pos(ind_z) >= corner(ind_z) .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z) / 2.D0) THEN
 !   print*, 'next_cell_up: cross 2'
     if(next_leak == posy) act_cell = upper_cell + 1
     if(next_leak == negy) act_cell = upper_cell + 3
    ! upper left cell
    ELSE IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) / 2.D0 .AND. &
       cross_pos(ind_z) >= corner(ind_z) + width(ind_z) / 2.D0 .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z)) THEN
 !   print*, 'next_cell_up: cross 3'
     if(next_leak == posy) act_cell = upper_cell + 4
     if(next_leak == negy) act_cell = upper_cell + 6
    ! upper right cell
    ELSE IF(cross_pos(ind_x) >= corner(ind_x) + width(ind_x) / 2.D0 .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) .AND. &
       cross_pos(ind_z) >= corner(ind_z) + width(ind_z) / 2.D0 .AND. cross_pos(ind_z) <= corner(ind_z) + width(ind_z)) THEN
 !   print*, 'next_cell_up: cross 4'
     if(next_leak == posy) act_cell = upper_cell + 5
     if(next_leak == negy) act_cell = upper_cell + 7
    ELSE
     write(*,*) 'next_cell_up: no cell was found'
     write(*,*) 'pos = ', norm2(cross_pos) / R_inf
     CALL abort()
    END IF
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! cross in z direction
   ELSE IF(next_leak == posz .OR. next_leak == negz) THEN
 !   print*, 'next_cell_up: cross = ', cross
    ! front left cell
    IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) / 2.D0 .AND. &
       cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y) / 2.D0) THEN
 !   print*, 'next_cell_up: cross 1'
     if(next_leak == posz) act_cell = upper_cell
     if(next_leak == negz) act_cell = upper_cell + 4
    ! lower right cell
    ELSE IF(cross_pos(ind_x) >= corner(ind_x) + width(ind_x) / 2.D0 .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) .AND. &
       cross_pos(ind_y) >= corner(ind_y) .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y) / 2.D0) THEN
 !   print*, 'next_cell_up: cross 2'
     if(next_leak == posz) act_cell = upper_cell + 1
     if(next_leak == negz) act_cell = upper_cell + 5
    ! upper left cell
    ELSE IF(cross_pos(ind_x) >= corner(ind_x) .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) / 2.D0 .AND. &
       cross_pos(ind_y) >= corner(ind_y) + width(ind_y) / 2.D0 .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y)) THEN
 !   print*, 'next_cell_up: cross 3'
     if(next_leak == posz) act_cell = upper_cell + 2
     if(next_leak == negz) act_cell = upper_cell + 6
    ! upper right cell
    ELSE IF(cross_pos(ind_x) >= corner(ind_x) + width(ind_x) / 2.D0 .AND. cross_pos(ind_x) <= corner(ind_x) + width(ind_x) .AND. &
       cross_pos(ind_y) >= corner(ind_y) + width(ind_y) / 2.D0 .AND. cross_pos(ind_y) <= corner(ind_y) + width(ind_y)) THEN
 !   print*, 'next_cell_up: cross 4'
     if(next_leak == posz) act_cell = upper_cell + 3
     if(next_leak == negz) act_cell = upper_cell + 7
    ELSE
     write(*,*) 'next_cell_up: no cell was found'
     write(*,*) 'pos = ', norm2(cross_pos)
     write(*,*) 'cor_z = ', corner(ind_z)/R_star, ' pos_z = ', norm2(cross_pos), ' cor_z+w = ', (corner(ind_z) + &
     width(ind_z))/R_star
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
  act_cell = next_cell
  upper_cell = dyn_cell(act_cell)%up_cell
  IF(upper_cell == 0) THEN
   next_cell = act_cell
   ! RETURN POINT
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
 
  IF(next_leak == posx) THEN
   subind_x = 1
  ELSE IF(next_leak == negx) THEN
   subind_x = sub_nx
  ELSE IF(next_leak == posy) THEN
   subind_y = 1
  ELSE IF(next_leak == negy) THEN
   subind_y = sub_ny
  ELSE IF(next_leak == posz) THEN
   subind_z = 1
  ELSE IF(next_leak == negz) THEN
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







END IF










END SUBROUTINE d_choosenextcell
