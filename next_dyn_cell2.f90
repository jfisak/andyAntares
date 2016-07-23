! this subroutine calculates the next dynamic cell
! which will be crossed by photon
  SUBROUTINE next_dyn_cell(pos,act_cell,pack_index,next_cell)

  USE types
  IMPLICIT NONE
  
  ! input variables
  INTEGER                               :: act_cell, pack_index
  ! output
  INTEGER                               :: next_cell
  DOUBLE PRECISION, DIMENSION(3)        :: pos
  ! number of bounds in the basic cell
  DOUBLE PRECISION, DIMENSION(3)        :: basic_width, ph_pos
  INTEGER                               :: bound_bc
  ! basic cell variables
  INTEGER                               :: ind_x, ind_y, ind_z, ind_cell_numb
  INTEGER                               :: bound1, bound2, bound3, bound4
  INTEGER                               :: bound5, bound6, bound7, bound8
  INTEGER                               :: bound21, bound22, bound23, bound24
  INTEGER                               :: nx_bound, ny_bound, nz_bound
  ! next basic cell
  INTEGER                               :: next_bcell

  basic_width = dyn_cell(1)%width
  ph_pos = package(pack_index)%pos
  bound_bc = 0
  ! firstly we have to compute the position in the basic cell
  IF(FLOOR((pos(1) + xmax)/basic_width(1)) == ((pos(1) + xmax))/basic_width(1)) THEN
   bound_bc = bound_bc + 1
   nx_bound = (pos(1) + xmax) / basic_width(1) + 1
  ELSE
   nx_bound = 0
  END IF
  IF(FLOOR((pos(2) + ymax)/basic_width(2)) == ((pos(2) + ymax))/basic_width(2)) THEN
   bound_bc = bound_bc + 1
   ny_bound = (pos(2) + ymax) / basic_width(2) + 1
  ELSE
   ny_bound = 0
  END IF
  IF(FLOOR((pos(3) + zmax)/basic_width(3)) == ((pos(3) + zmax))/basic_width(3)) THEN
   bound_bc = bound_bc + 1
   nz_bound = (pos(3) + zmax) / basic_width(3) + 1
  ELSE
   nz_bound = 0
  END IF

  ! we have to find all basic cells containing this point
  ! bound_bc = 0 inside the basic cell => we have to look only into one bc
  ! bound_bc = 1 in the wall of bc
  ! bound_bc = 2 in the edge of bc
  ! bound_bc = 3 in the corner of bc
  print*, 'bound_bc = ', bound_bc
  ! firstly we have to know the basic cell of the photon
     ind_x = FLOOR(ph_pos(1)/basic_width(1) + DBLE(nx_cell)/2) + 1
     ind_y = FLOOR(ph_pos(2)/basic_width(2) + DBLE(ny_cell)/2) + 1
     ind_z = FLOOR(ph_pos(3)/basic_width(3) + DBLE(nz_cell)/2) + 1
     ind_cell_numb = (ind_x - 1) * ny_cell * nz_cell + (ind_y - 1) * nz_cell + ind_z
  ! move in the same basic cell
  IF(bound_bc == 0) THEN
   next_bcell = ind_cell_numb
  ! moving to another basic cell
  ELSE IF(bound_bc == 1) THEN
   ! will the photon stay in the grid?
   IF(nx_bound == 1 .OR. nx_bound == xmax + 1 .OR. &
      ny_bound == 1 .OR. ny_bound == ymax + 1 .OR. &
      nz_bound == 1 .OR. nz_bound == zmax + 1) THEN
      next_cell = -99
   END IF
   ! if it stays in grid we have to calculate
   ! index of next cell
   IF(nx_bound /= 0) THEN
    bound1 = (nx_bound - 2) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
    if(bound1 /= ind_cell_numb) next_cell = bound1
    if(bound2 /= ind_cell_numb) next_cell = bound2
   ELSE IF(ny_bound /= 0) THEN
    bound1 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 2) * nz_cell + nz_bound
    bound2 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
    if(bound1 /= ind_cell_numb) next_cell = bound1
    if(bound2 /= ind_cell_numb) next_cell = bound2
   ELSE IF(nz_bound /= 0) THEN
    bound1 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + (nz_bound - 1)
    bound2 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
    if(bound1 /= ind_cell_numb) next_bcell = bound1
    if(bound2 /= ind_cell_numb) next_bcell = bound2
   END IF 
   ! moving to another basic cell through an edge
  ELSE IF(bound_bc == 2) THEN
   IF(nx_bound /= 0) THEN
    bound1 = (nx_bound - 2) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
    bound2 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
   ELSE IF(ny_bound /= 0 .AND. nx_bound == 0) THEN
    bound1 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 2) * nz_cell + nz_bound
    bound2 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
   ELSE IF(ny_bound /= 0 .AND. nx_bound /= 0) THEN
    bound3 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 2) * nz_cell + nz_bound
    bound4 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
   ELSE IF(nz_bound /= 0) THEN
    bound3 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + (nz_bound - 1)
    bound4 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
   END IF 
   ! now we have to find which two indeces are same
   ! because bound4 = bound2
   bound21 = bound4
   bound22 = bound1 + bound3 - bound2
   bound23 = bound1
   bound24 = bound3
   IF(bound21 == ind_cell_numb) next_bcell = bound24
   IF(bound22 == ind_cell_numb) next_bcell = bound22
   IF(bound23 == ind_cell_numb) next_bcell = bound23
   IF(bound24 == ind_cell_numb) next_bcell = bound21

   ! because ind_cell_numb + next_cell = bound22 + bound23
   next_cell = bound21 + bound22 - ind_cell_numb
  ELSE IF(bound_bc == 3) THEN
   bound1 = (nx_bound - 2) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound
   bound2 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 2) * nz_cell + nz_bound
   bound3 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + (nz_bound - 1)
   ! common cell
   bound4 = (nx_bound - 1) * ny_cell * nz_cell + (ny_bound - 1) * nz_cell + nz_bound

   bound5 = bound1 + bound2 - bound4
   bound6 = bound1 + bound3 - bound4
   bound7 = bound2 + bound3 - bound4
   bound8 = bound5 + bound6 - bound1
   
   IF(bound1 == ind_cell_numb) next_bcell = bound7
   IF(bound2 == ind_cell_numb) next_bcell = bound6
   IF(bound3 == ind_cell_numb) next_bcell = bound5
   IF(bound4 == ind_cell_numb) next_bcell = bound8
   IF(bound5 == ind_cell_numb) next_bcell = bound3
   IF(bound6 == ind_cell_numb) next_bcell = bound2
   IF(bound7 == ind_cell_numb) next_bcell = bound1
   IF(bound8 == ind_cell_numb) next_cell = bound4
  END IF

  
   

   

  
  


  END SUBROUTINE
