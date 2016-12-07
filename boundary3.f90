SUBROUTINE boundary3(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

  USE types

  IMPLICIT NONE    

INTEGER                         :: pack_index, next_cell
DOUBLE PRECISION                :: dist
! local actual cell
INTEGER                         :: actCell
! basic cell variables
! basic_cell_numb: actual basic cell number
INTEGER                         :: basic_cell_numb, next_bcell, next_bas_cell
DOUBLE PRECISION                :: bdist
! dynamic cell variables
INTEGER                         :: cell_numb
! photon properties
DOUBLE PRECISION, DIMENSION(3)  :: dir, phot_pos
! position of the cross point
DOUBLE PRECISION, DIMENSION(3)  :: cross_pos

! calculation of a distance from the basic cell
! firstly we have to know which basic cell photon occupies
!   Number of the current cell
  cell_numb = package(pack_index)%cell_numb
! define a local variable actCell
  actCell = cell_numb
!  print*, 'boundary3: pack index = ', pack_index, ' actCell = ', actCell
! number of the corresponding basic cell
  DO WHILE(dyn_cell(actCell)%down_cell /= 0) 
   actCell = dyn_cell(actCell)%down_cell
  END DO
! now we know a basic cell number
  basic_cell_numb = actCell
!  print*, 'boundary3: basic_cell_numb = ', basic_cell_numb
! now we compute the nearest distance from the basic cell
 CALL find_bdist(pack_index, basic_cell_numb, bdist)
! now we are computing the nearest distance to the actuall dynamic cell
 CALL find_dist(pack_index,cell_numb, dist)
! now we decide if the basic cell number will change or not
!  print*, 'boundary3: bdist = ', bdist, ' dist = ', dist
  IF(bdist == dist) THEN
   ! we have to find a new cell
   CALL find_basic_cell(basic_cell_numb, pack_index,bdist, next_bcell)
   next_bas_cell = next_bcell
  ELSE
   next_bas_cell = basic_cell_numb
  END IF
! position of the point
 cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist

! now we calculate next dynamical cell
  !print*, 'boundary3: cross_pos = ', cross_pos, ' cell_numb = ', cell_numb, &
  !              ' next_bas_cell = ', next_bas_cell
  IF(next_bas_cell <= 0) THEN
   next_cell = -99
  ELSE
   CALL find_dyn_cell2(cross_pos, cell_numb, next_bas_cell, next_cell)
  END IF
  !print*, 'next_cell = ', next_cell
END SUBROUTINE boundary3
