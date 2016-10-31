SUBROUTINE boundary3(pack_index, dist, next_cell) 
! write nx_cell, ny_cell, nz_cell as global variable

! Calculate the shortest distance to the cell surface which photon will cross and return this 
! distance (dist), the cell the photon will go to (next_cell, -99 in case of photon leaving 
! the simulation grid(volume)), and make the photon remember which surface will cross (update 
! package(pack_index)%last_cross).

  USE types

  IMPLICIT NONE    

    INTEGER                         :: pack_index, nc, next_cell,forbidden, hit_surface
    DOUBLE PRECISION                :: dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz
    DOUBLE PRECISION                :: dist_minz, dist, surface_pos 

    hit_surface = 0

! calculation of a distance from the basic cell
! firstly we have to know which basic cell photon occupies
!   Number of the current cell
  nc = package(pack_index)%cell_numb
! define a local variable actCell
  actCell = nc
! number of the corresponding basic cell
  DO WHILE(dyn_cell(actCell)%down_cell /= 0) 
   actCell = dyn_cell(actCell)%down_cell
  END DO
! now we know a basic cell number
  basic_cell_number = actCell
! now we compute the nearest distance from the basic cell
  bcorner = dyn_cell(cell_numb)%corner
  bwidth = dyn_cell(cell_numb)%width
  phot_pos = package(pack_index)%pos
  dir = package(pack_index)%dir
 ! we will calculate parameters t1,...,t6
 IF(dir(1) /= 0) THEN
  bt1 = (dyn_cell(basic_cell_numb)%corner(1) - package(pack_index)%pos(1))/(package(pack_index)%dir(1))
  bt4 = (dyn_cell(basic_cell_numb)%corner(1) + dyn_cell(basic_cell_numb)%width(1) - package(pack_index)%pos(1))/(package(pack_index)%dir(1))
 ELSE
  bt1 = 0
  bt4 = 0
 END IF
 IF(dir(2) /= 0) THEN
  bt2 = (dyn_cell(basic_cell_numb)%corner(2) - package(pack_index)%pos(2))/(package(pack_index)%dir(2))
  bt5 = (dyn_cell(basic_cell_numb)%corner(2) + dyn_cell(basic_cell_numb)%width(2) - package(pack_index)%pos(2))/(package(pack_index)%dir(2))
 ELSE
  bt2 = 0
  bt5 = 0
 END IF
 IF(dir(3) /= 0) THEN
  bt3 = (dyn_cell(basic_cell_numb)%corner(3) - package(pack_index)%pos(3))/(package(pack_index)%dir(3))
  bt6 = (dyn_cell(basic_cell_numb)%corner(3) + dyn_cell(basic_cell_numb)%width(3) - package(pack_index)%pos(3))/(package(pack_index)%dir(3))
 ELSE
  bt3 = 0
  bt6 = 0
 END IF

  bdist = 1.D99
  ! we are looking for a bound in front of the photon,
  ! so we have to choose solution with t > 0
  IF( (bt1 > 0.E0) .AND. (bt1 < dist) ) THEN
   bdist = bt1
  END IF
  IF( (bt2 > 0.E0)  .AND. (bt2 < dist) ) THEN
   bdist = bt2
  END IF
  IF( (bt3 > 0.E0) .AND. (bt3 < dist)  ) THEN
   bdist = bt3
  END IF
  IF( (bt4 > 0.E0) .AND. (bt4 < dist)  ) THEN
   bdist = bt4
  END IF
  IF( (bt5 > 0.E0) .AND. (bt5 < dist)  ) THEN
   bdist = bt5
  END IF
  IF( (bt6 > 0.E0) .AND. (bt6 < dist)  ) THEN
   bdist = bt6
  END IF
  
 ! now we are computing the nearest distance to the actuall dynamic cell
  corner = dyn_cell(cell_numb)%corner
  width = dyn_cell(cell_numb)%width
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
    CALL next_dyn_cell3(nc,pack_index,next_cell)
    !print*, 'nc = ', nc, 'pack_index = ', pack_index, 'next_cell = ', next_cell, 'dist = ', dist
!   Calculate the distances to the all cell surfaces from the photon current position along the ray
!   (formula for this can be found in http://www.roe.ac.uk/ifa/postgrad/pedagogy/2009_forgan.pdf 
!    - Fig. 2, An Introduction to Monte Carlo Radiative Transfer, Duncan Forgan)
    IF (debug .EQ. 1) THEN 
        print*, 'boundary: dir ', package(pack_index)%dir
        print*, 'boundary: distances posx, negx, posy, negy, posz, negz ', &
             dist_plusx, dist_minx, dist_plusy, dist_miny, dist_plusz, dist_minz
    END IF
   
END SUBROUTINE boundary3
