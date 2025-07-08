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
DOUBLE PRECISION, DIMENSION(3)          :: cross_pos
INTEGER                                 :: cross
! actual cell
INTEGER                                 :: act_cell
DOUBLE PRECISION, DIMENSION(3)          :: corner, upcorner, centre
INTEGER                                 :: upper_cell
! parameters of subcells of dyngrid ijk
DOUBLE PRECISION, DIMENSION(3)          :: subcells_width
INTEGER                                 :: subind_x, subind_y, subind_z
INTEGER                                         :: sub_nx, sub_ny, sub_nz

INTEGER, DIMENSION(3)                   :: bcell
INTEGER                                 :: bindex
DOUBLE PRECISION, DIMENSION(3)          :: pos

DOUBLE PRECISION, PARAMETER                     :: epsilon0 = 1e-6
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
 centre = (upcorner - corner) / 2.D0
 upcorner = dyn_cell(act_cell)%upcorner
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
    IF(cross_pos(2) >= upcorner(2) .AND. cross_pos(2) <= centre(2)  .AND. &
       cross_pos(3) >= upcorner(3) .AND. cross_pos(3) <= centre(3) ) THEN
 !   print*, 'next_cell_up: cross 1'
     if(cross == posx) act_cell = upper_cell
     if(cross == negx) act_cell = upper_cell + 1
    ! lower rear cell
    ELSE IF(cross_pos(2) >= centre(2) .AND. cross_pos(2) <= upcorner(2) .AND. &
            cross_pos(3) >= corner(3) .AND. cross_pos(3) <= centre(3)) THEN
 !   print*, 'next_cell_up: cross 2'
     if(cross == posx) act_cell = upper_cell + 2
     if(cross == negx) act_cell = upper_cell + 3
    ! upper front cell
    ELSE IF(cross_pos(2) >= corner(2) .AND. cross_pos(2) <= centre(2) .AND. &
       cross_pos(3) >= centre(3) .AND. cross_pos(3) <= upcorner(3)) THEN
 !   print*, 'next_cell_up: cross 3'
     if(cross == posx) act_cell = upper_cell + 4
     if(cross == negx) act_cell = upper_cell + 5
    ! upper rear cell
    ELSE IF(cross_pos(2) >= centre(2) .AND. cross_pos(2) <= upcorner(2) .AND. &
       cross_pos(3) >= centre(3) .AND. cross_pos(3) <= upcorner(3)) THEN
 !   print*, 'next_cell_up: cross 4'
     if(cross == posx) act_cell = upper_cell + 6
     if(cross == negx) act_cell = upper_cell + 7
    ELSE
     ! write(*,*) 'next_cell_up: no cell was found'
     ! CALL abort()
    END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! cross in y direction
  ELSE IF(cross == posy .OR. cross == negy) THEN
!   print*, 'next_cell_up: cross = ', cross
   ! lower left cell
   IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= centre(1)  .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= centre(3) ) THEN
!   print*, 'next_cell_up: cross 1'
    if(cross == posy) act_cell = upper_cell
    if(cross == negy) act_cell = upper_cell + 2
   ! lower right cell
   ELSE IF(cross_pos(1) >= centre(1) .AND. cross_pos(1) <= upcorner(1) .AND. &
      cross_pos(3) >= corner(3) .AND. cross_pos(3) <= centre(3)) THEN
!   print*, 'next_cell_up: cross 2'
    if(cross == posy) act_cell = upper_cell + 1
    if(cross == negy) act_cell = upper_cell + 3
   ! upper left cell
   ELSE IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= centre(1) .AND. &
      cross_pos(3) >= centre(3) .AND. cross_pos(3) <= upcorner(3)) THEN
!   print*, 'next_cell_up: cross 3'
    if(cross == posy) act_cell = upper_cell + 4
    if(cross == negy) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(1) >= upcorner(1) / 2.D0 .AND. cross_pos(1) <= upcorner(1) .AND. &
      cross_pos(3) >= centre(3) .AND. cross_pos(3) <= upcorner(3)) THEN
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
   IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= centre(1) .AND. &
      cross_pos(2) >= corner(2) .AND. cross_pos(2) <= centre(2)) THEN
!   print*, 'next_cell_up: cross 1'
    if(cross == posz) act_cell = upper_cell
    if(cross == negz) act_cell = upper_cell + 4
   ! lower right cell
   ELSE IF(cross_pos(1) >= centre(1) .AND. cross_pos(1) <= upcorner(1) .AND. &
      cross_pos(2) >= corner(2) .AND. cross_pos(2) <= centre(2)) THEN
!   print*, 'next_cell_up: cross 2'
    if(cross == posz) act_cell = upper_cell + 1
    if(cross == negz) act_cell = upper_cell + 5
   ! upper left cell
   ELSE IF(cross_pos(1) >= corner(1) .AND. cross_pos(1) <= centre(1) .AND. &
      cross_pos(2) >= centre(2) .AND. cross_pos(2) <= centre(2)) THEN
!   print*, 'next_cell_up: cross 3'
    if(cross == posz) act_cell = upper_cell + 2
    if(cross == negz) act_cell = upper_cell + 6
   ! upper right cell
   ELSE IF(cross_pos(1) >= centre(1) .AND. cross_pos(1) <= centre(1) .AND. &
      cross_pos(2) >= centre(2) .AND. cross_pos(2) <= centre(2)) THEN
!   print*, 'next_cell_up: cross 4'
    if(cross == posz) act_cell = upper_cell + 3
    if(cross == negz) act_cell = upper_cell + 7
   ELSE
    write(*,*) 'next_cell_up: no cell was found'
    write(*,*) 'pos = ', norm2(package(pack_index)%pos) / R_star
    ! write(*,*) 'cor_z = ', corner(3)/R_star, ' pos_z = ', package(pack_index)%pos(3)/R_star, ' cor_z+w = ', (corner(3) + &
    ! width(3))/R_star
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
 subcells_width = dyn_cell(upper_cell)%width
 subcells_width = dyn_cell(upper_cell)%width
 sub_nx = dyn_cell(act_cell)%n_sbgr(1)
 sub_ny = dyn_cell(act_cell)%n_sbgr(2)
 sub_nz = dyn_cell(act_cell)%n_sbgr(3)

 subind_x = FLOOR((cross_pos(1) - dyn_cell(act_cell)%corner(1))/subcells_width(1)) + 1
 subind_y = FLOOR((cross_pos(2) - dyn_cell(act_cell)%corner(2))/subcells_width(2)) + 1
 subind_z = FLOOR((cross_pos(3) - dyn_cell(act_cell)%corner(3))/subcells_width(3)) + 1

 pos = package(pack_index)%pos
 bcell(1) = FLOOR(pos(1)/basic_cell_width(1) + dble(nx_cell)/2.D0) + 1
 bcell(2) = FLOOR(pos(2)/basic_cell_width(2) + dble(ny_cell)/2.D0) + 1
 bcell(3) = FLOOR(pos(3)/basic_cell_width(3) + dble(nz_cell)/2.D0) + 1
  
 ! index of the given basic cell
 bindex = (bcell(1) - 1) * ny_cell * nz_cell + (bcell(2) - 1) * nz_cell + bcell(3)
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
