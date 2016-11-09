SUBROUTINE boundary3(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

  USE types

  IMPLICIT NONE    

    INTEGER                         :: pack_index, next_cell,forbidden, hit_surface
    DOUBLE PRECISION                :: dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz
    DOUBLE PRECISION                :: dist_minz, dist, surface_pos 
    ! photon properties
    DOUBLE PRECISION, DIMENSION(3)  :: dir, phot_pos
    ! variables for dynamic cells
    INTEGER                         :: cell_numb, act_dyn_cell
    DOUBLE PRECISION                :: t1, t2, t3, t4, t5, t6
    DOUBLE PRECISION, DIMENSION(3)  :: corner, width
    ! variables for basic cells
    INTEGER                         :: actCell, basic_cell_numb
    INTEGER                         :: next_bas_cell, next_bcell
    DOUBLE PRECISION, DIMENSION(3)  :: bcorner, bwidth
    DOUBLE PRECISION                :: bt1, bt2, bt3, bt4, bt5, bt6, bdist
    ! position of the cross point
    DOUBLE PRECISION, DIMENSION(3)  :: cross_pos

    hit_surface = 0

! calculation of a distance from the basic cell
! firstly we have to know which basic cell photon occupies
!   Number of the current cell
  cell_numb = package(pack_index)%cell_numb
! define a local variable actCell
  actCell = cell_numb
  print*, 'boundary3: actCell = ', actCell
! number of the corresponding basic cell
  DO WHILE(dyn_cell(actCell)%down_cell /= 0) 
   actCell = dyn_cell(actCell)%down_cell
  END DO
! now we know a basic cell number
  basic_cell_numb = actCell
! now we compute the nearest distance from the basic cell
  bcorner = dyn_cell(cell_numb)%corner
  bwidth = dyn_cell(cell_numb)%width
  phot_pos = package(pack_index)%pos
  dir = package(pack_index)%dir
 ! we will calculate parameters t1,...,t6
 IF(dir(1) /= 0) THEN
  bt1 = (dyn_cell(basic_cell_numb)%corner(1) - package(pack_index)%pos(1))/(package(pack_index)%dir(1))
  bt4 = (dyn_cell(basic_cell_numb)%corner(1) + dyn_cell(basic_cell_numb)%width(1) - &
        package(pack_index)%pos(1))/(package(pack_index)%dir(1))
 ELSE
  bt1 = 0
  bt4 = 0
 END IF
 IF(dir(2) /= 0) THEN
  bt2 = (dyn_cell(basic_cell_numb)%corner(2) - package(pack_index)%pos(2))/(package(pack_index)%dir(2))
  bt5 = (dyn_cell(basic_cell_numb)%corner(2) + dyn_cell(basic_cell_numb)%width(2) - &
        package(pack_index)%pos(2))/(package(pack_index)%dir(2))
 ELSE
  bt2 = 0
  bt5 = 0
 END IF
 IF(dir(3) /= 0) THEN
  bt3 = (dyn_cell(basic_cell_numb)%corner(3) - package(pack_index)%pos(3))/(package(pack_index)%dir(3))
  bt6 = (dyn_cell(basic_cell_numb)%corner(3) + dyn_cell(basic_cell_numb)%width(3) - &
        package(pack_index)%pos(3))/(package(pack_index)%dir(3))
 ELSE
  bt3 = 0
  bt6 = 0
 END IF

  bdist = 1.D99
  ! we are looking for a bound in front of the photon,
  ! so we have to choose solution with t > 0
  IF( (bt1 > 0.E0) .AND. (bt1 < dist) ) THEN
   bdist = bt1
   IF(dyn_cell(basic_cell_numb)%indexc(1) == nx_cell) THEN
    next_bcell = -99
   ELSE
    next_bcell = basic_cell_numb + ny_cell * nz_cell
   END IF
  END IF
  IF( (bt2 > 0.E0) .AND. (bt2 < dist) ) THEN
   bdist = bt2
   IF(dyn_cell(basic_cell_numb)%indexc(1) == ny_cell) THEN
    next_bcell = -99
   ELSE
    next_bcell = basic_cell_numb + nz_cell
   END IF
  END IF
  IF( (bt3 > 0.E0) .AND. (bt3 < dist)  ) THEN
   bdist = bt3
   IF(dyn_cell(basic_cell_numb)%indexc(1) == nz_cell) THEN
    next_bcell = -99
   ELSE
    next_bcell = basic_cell_numb + ny_cell * nz_cell
   END IF
  END IF
  IF( (bt4 > 0.E0) .AND. (bt4 < dist)  ) THEN
   bdist = bt4
   IF(dyn_cell(basic_cell_numb)%indexc(1) == 1) THEN
    next_bcell = -99
   ELSE
    next_bcell = basic_cell_numb - ny_cell * nz_cell
   END IF
  END IF
  IF( (bt5 > 0.E0) .AND. (bt5 < dist)  ) THEN
   bdist = bt5
   IF(dyn_cell(basic_cell_numb)%indexc(1) == 1) THEN
    next_bcell = -99
   ELSE
    next_bcell = basic_cell_numb - nz_cell
   END IF
  END IF
  IF( (bt6 > 0.E0) .AND. (bt6 < dist)  ) THEN
   bdist = bt6
   IF(dyn_cell(basic_cell_numb)%indexc(1) == 1) THEN
    next_bcell = -99
   ELSE
    next_bcell = basic_cell_numb - 1
   END IF
  END IF
  
 ! now we are computing the nearest distance to the actuall dynamic cell
  corner = dyn_cell(cell_numb)%corner
  width = dyn_cell(cell_numb)%width
  act_dyn_cell = package(pack_index)%cell_numb
 ! we will calculate parameters t1,...,t6
 IF(dir(1) /= 0) THEN
  t1 = (dyn_cell(cell_numb)%corner(1) - package(pack_index)%pos(1))/(package(pack_index)%dir(1))
  t4 = (dyn_cell(cell_numb)%corner(1) + dyn_cell(cell_numb)%width(1) - package(pack_index)%pos(1))/(package(pack_index)%dir(1))
 ELSE
  t1 = 0
  t4 = 0
 END IF
 IF(dir(2) /= 0) THEN
  t2 = (dyn_cell(cell_numb)%corner(2) - package(pack_index)%pos(2))/(package(pack_index)%dir(2))
  t5 = (dyn_cell(cell_numb)%corner(2) + dyn_cell(cell_numb)%width(2) - package(pack_index)%pos(2))/(package(pack_index)%dir(2))
 ELSE
  t2 = 0
  t5 = 0
 END IF
 IF(dir(3) /= 0) THEN
  t3 = (dyn_cell(cell_numb)%corner(3) - package(pack_index)%pos(3))/(package(pack_index)%dir(3))
  t6 = (dyn_cell(cell_numb)%corner(3) + dyn_cell(cell_numb)%width(3) - package(pack_index)%pos(3))/(package(pack_index)%dir(3))
 ELSE
  t3 = 0
  t6 = 0
 END IF

  dist = 1.D99
  ! we are looking for the bound in front of the photon,
  ! so we have to choose solution with t > 0
  IF( (t1 > 0.E0) .AND. (t1 < dist) ) THEN
   dist = t1
  END IF
  IF( (t2 > 0.E0)  .AND. (t2 < dist) ) THEN
   dist = t2
  END IF
  IF( (t3 > 0.E0) .AND. (t3 < dist)  ) THEN
   dist = t3
  END IF
  IF( (t4 > 0.E0) .AND. (t4 < dist)  ) THEN
   dist = t4
  END IF
  IF( (t5 > 0.E0) .AND. (t5 < dist)  ) THEN
   dist = t5
  END IF
  IF( (t6 > 0.E0) .AND. (t6 < dist)  ) THEN
   dist = t6
  END IF
  ! now we decide if the basic cell number will change or not
  IF(bdist == dist) THEN
   next_bas_cell = next_bcell
  ELSE
   next_bas_cell = basic_cell_numb
  END IF

  ! position of the point
  cross_pos = package(pack_index)%pos + package(pack_index)%dir * dist

  ! now we calculate next dynamical cell
  IF(next_bas_cell /= -99) CALL find_dyn_cell2(cross_pos, act_dyn_cell, next_bas_cell, next_cell)
   
END SUBROUTINE boundary3
