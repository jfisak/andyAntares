! this subroutine will create a dynamical grid cell for the given basic cell by
! calling subroutine divide_cell_8 or divide_cell_ijk
!
! INPUT: n_dyncell(INT): index of the upper propGrid cell
!        max_n_dcell(INT): the total number of created propGrid cells
!
  SUBROUTINE create_dynamical_grid_cells(n_dyncell,max_n_dcell)

   USE types
USE constants

   IMPLICIT NONE
   
   INTEGER                              :: n_dyncell
   ! TYPE(virt_point), DIMENSION(Npart):: v_part
   ! number of created dynamic cells and
   ! actual number of grid cell
   INTEGER                              :: max_n_dcell, act_n_dyncell
   INTEGER                              :: ind_I, ind_J
   INTEGER                              :: n_points
   ! number of new created cells in cell
   INTEGER                              :: no_dcells
   ! variables for boundaries
   INTEGER                              :: up_bound, newbound
   ! maximal number of point in one cell
   INTEGER, PARAMETER                   :: maxPart = 1
   DOUBLE PRECISION, DIMENSION(const_dimofspace)       :: corner, cell_width_2
   TYPE(virt_point), ALLOCATABLE     :: local_point(:)
   TYPE(dyn_grid_cell), ALLOCATABLE     ::  pom2(:)
   ! dimension of the subcell grid
   INTEGER, DIMENSION(const_dimofspace)                :: dimofsubcells
   ! move to the next cell, if the size is smaller than minimal possible cell size
   LOGICAL                              :: next_cell
   ! local variables
   DOUBLE PRECISION, DIMENSION(const_dimofspace)       :: loc_corner, loc_cell_width
   DOUBLE PRECISION, DIMENSION(const_dimofspace)       :: vp_pos
   INTEGER                              :: loc_np
   INTEGER                              :: loc_downcell, loc_upcell
   INTEGER                              :: cur_point
    ! for 8-dyncells

corner(:) = dyn_cell(n_dyncell)%corner(:)
cell_width_2(:) = dyn_cell(n_dyncell)%width(:)


! at first we have to know, how many points are located
! in the given cell
n_points = dyn_cell(n_dyncell)%n_virt
SELECT CASE(dyngrid)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
! 8-type dynamical cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
CASE(1)
ALLOCATE(local_point(n_points))

cur_point = 0
DO ind_I = 1, Nvirtpoint
 vp_pos = virtual_point(ind_I)%pos
 IF(vp_pos(ind_x) > corner(ind_x) .and. vp_pos(ind_x) < (corner(ind_x) + cell_width_2(ind_x)) .and. &
  vp_pos(ind_y) >= corner(ind_y) .and. vp_pos(ind_y) < (corner(ind_y) + cell_width_2(ind_y)) .and. &
  vp_pos(ind_z) >= corner(ind_z) .and. vp_pos(ind_z) < (corner(ind_z) + cell_width_2(ind_z))) THEN
  cur_point = cur_point + 1
  local_point(cur_point) = virtual_point(ind_I)
 END IF
END DO

! we have found point included in this basic cell
! now we have to generate brand new dynamical cell
! if the number of point n_points is equal to one or two
! we don't have to allocate any new dynamical cells
! number of the dynamic grid cell
! firstly for the basic cell division
!n_dg = 0
up_bound = SIZE(dyn_cell(:))
act_n_dyncell = n_dyncell
no_dcells = 8
! start: large loop
DO
! print*, 'create_dynamical_grid_cells: act_n_dyncell = ', act_n_dyncell
! local number of point is in the begining of cycle = 0
  loc_np = 0
  loc_corner(:) = dyn_cell(act_n_dyncell)%corner(:)
  loc_cell_width(:) = dyn_cell(act_n_dyncell)%width(:)
  loc_upcell = dyn_cell(act_n_dyncell)%up_cell
  loc_downcell = dyn_cell(act_n_dyncell)%down_cell
 ! how many virtual point is there in this subcell
 IF(next_cell .EQV. .FALSE.) THEN
  ! start: calculating number of local point
  DO ind_J = 1, n_points
   IF((local_point(ind_J)%pos(ind_x) >= loc_corner(ind_x)) .AND. &
     (local_point(ind_J)%pos(ind_x) < (loc_corner(ind_x) + loc_cell_width(ind_x))) .AND. &
     (local_point(ind_J)%pos(ind_y) >= loc_corner(ind_y)) .AND. &
     (local_point(ind_J)%pos(ind_y) < (loc_corner(ind_y) + loc_cell_width(ind_y))) .AND. &
     (local_point(ind_J)%pos(ind_z) >= loc_corner(ind_z)) .AND. &
     (local_point(ind_J)%pos(ind_z) < (loc_corner(ind_z) + loc_cell_width(ind_z)))) THEN
   !print*, 'found a point number ', loc_n_points + 1
   loc_np = loc_np + 1
  END IF
 ! stop: calculating number of local point
  END DO
 ELSE
  loc_np = 0
  next_cell = .FALSE.
 END IF
 ! now we have to decide what to do on the basement of number of
 ! local point
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
 ! 1. the division of the cells is good enough
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!§§
  !print*, 'dynamic grid: ', up_bound, act_n_dyncell, max_n_dcell
 IF(loc_np <= maxPart ) THEN
  ! start: we move to the lower level of the grid to the "not ending" subcell
  DO
!    IF(loc_downcell /= 0) &
!      print*, 'dyn_cell(down_cell)%up_cell - act_n_dyncell == no_dcells - 1', &
!      act_n_dyncell - dyn_cell(loc_downcell)%up_cell 
   IF(dyn_cell(act_n_dyncell)%down_cell /= 0) THEN
     IF(act_n_dyncell - dyn_cell(loc_downcell)%up_cell == no_dcells - 1) &
      THEN
      ! we will move one level lower
      !print*, 'moving to the lower level...'
      act_n_dyncell = dyn_cell(act_n_dyncell)%down_cell
      loc_upcell = dyn_cell(act_n_dyncell)%up_cell
      loc_downcell = dyn_cell(act_n_dyncell)%down_cell
     ! we can go to the next subcell in the given level
     ELSE
     ! print*, 'moving to the next cell...'
      act_n_dyncell = act_n_dyncell + 1
      EXIT
     END IF
   ELSE
    EXIT
   END IF
  ! stop: we move to the lower level of the grid to the "not ending" subcell
  END DO
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 2. we have to create additional cells otherwise
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ELSE
  ! we check if the width of new cells will be large enough
  if (dyn_cell(act_n_dyncell)%width(ind_x)/2.E0 < minwidth .OR. &
   dyn_cell(act_n_dyncell)%width(ind_y)/2.E0 < minwidth .OR. &
       dyn_cell(act_n_dyncell)%width(ind_z)/2.E0 < minwidth) then
       loc_np = 0
      ! print*, 'CELL WOULD BE TOO SMALL...MOVING TO THE NEXT CELL...'
       next_cell = .TRUE.
   end if
  ! is there some free space left in the field dyn_cell?
   up_bound = SIZE(dyn_cell(:))
   do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do we need to resize dyn_cell?
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if(max_n_dcell + no_dcells >= up_bound) then
     ! define a new upper bound
     newbound = 2 * up_bound
     ALLOCATE(pom2(up_bound))
     do ind_I = 1, up_bound
      pom2(ind_I) = dyn_cell(ind_I)
     end do
     DEALLOCATE(dyn_cell)
     ALLOCATE(dyn_cell(newbound))
     do ind_I = 1, up_bound
      dyn_cell(ind_I) = pom2(ind_I)
     end do
     DEALLOCATE(pom2)
     up_bound = newbound
    else
     exit
    end if
   end do
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! now we will create new subcells
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   CALL divide_cell_8(act_n_dyncell, max_n_dcell)
   ! we have to increase the variable max_n_dcell
   max_n_dcell = max_n_dcell + no_dcells
   act_n_dyncell = dyn_cell(act_n_dyncell)%up_cell
   ! if every single cell is divided correctly the loop will stop
  END IF
 ! can we stop the large loop?
 ! only if we come back to the basic cell we are creating subcells in
 IF(act_n_dyncell == n_dyncell) THEN
!  print*, 'all dynamic cells for this basic cell were created...'
  EXIT
 END IF
! stop: large loop
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ijk-type dynamical cell
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE(2)
 IF(model_type == 1) THEN
  dimofsubcells(ind_x) = FLOOR(n_points**(0.5))
  dimofsubcells(ind_y) = FLOOR(n_points**(0.5))
  dimofsubcells(ind_z) = FLOOR(n_points**(0.5))
 ELSE IF(model_type == 3) THEN
  dimofsubcells(ind_x) = FLOOR(n_points**(0.5))
  dimofsubcells(ind_y) = FLOOR(n_points**(0.5))
  dimofsubcells(ind_z) = FLOOR(n_points**(0.5))
 ELSE
  dimofsubcells(ind_x) = FLOOR(n_points**(4.0))
  dimofsubcells(ind_y) = FLOOR(n_points**(4.0))
  dimofsubcells(ind_z) = FLOOR(n_points**(4.0))
 END IF
  dyn_cell(n_dyncell)%n_sbgr = dimofsubcells
  no_dcells = dimofsubcells(ind_x) * dimofsubcells(ind_y) * dimofsubcells(ind_z)
 ! write(*,*) 'create_dynamical_grid_cells: no_dcells = ', no_dcells
 IF(dimofsubcells(ind_x) <= 1 .and. dimofsubcells(ind_y) <= 1 &
  .and. dimofsubcells(ind_z) <=1) THEN
  no_dcells = 0
 END IF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! do we need to resize dyn_cell?
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  up_bound = SIZE(dyn_cell(:))
  ! write(*,*) 'create_dynamical_grid_cells: max_n_dcell = ', max_n_dcell, ' up_bound = ', up_bound
  if(max_n_dcell + no_dcells >= up_bound) then
  ! write(*,*)  'create_dynamical_grid_cells: creating a larger array dyn_cell...'
   ! define a new upper bound
   newbound =  max_n_dcell + 2*no_dcells
   ALLOCATE(pom2(up_bound))
   do ind_I = 1, up_bound
    pom2(ind_I) = dyn_cell(ind_I)
   end do
   DEALLOCATE(dyn_cell)
   ALLOCATE(dyn_cell(newbound))
   do ind_I = 1, up_bound
    dyn_cell(ind_I) = pom2(ind_I)
   end do
   DEALLOCATE(pom2)
   up_bound = newbound
  end if
  IF(n_points >= 8 .OR. no_dcells > 1) THEN
   CALL divide_cell_ijk(n_dyncell, max_n_dcell, dimofsubcells)
   max_n_dcell = max_n_dcell + no_dcells
  END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! default case
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CASE DEFAULT
 STOP 'create_dynamical_grid_cells: choice of the dyngrid type is not known'
END SELECT

END SUBROUTINE create_dynamical_grid_cells
