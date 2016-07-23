! this subroutine find
  SUBROUTINE next_dyn_cell3(cell_numb,pack_index,next_cell,dist)

  USE types
  IMPLICIT NONE

  ! input variables
  INTEGER                               :: cell_numb, pack_index
  ! output variables
  DOUBLE PRECISION, DIMENSION(3)        :: pos
  INTEGER                               :: next_cell
  ! parameter
  DOUBLE PRECISION                      :: t1, t2, t3, t4, t5, t6
  DOUBLE PRECISION                      :: tTest, tbound
  DOUBLE PRECISION                      :: solution, dist
  ! 
  DOUBLE PRECISION, DIMENSION(3)        :: corner, width, phot_pos, dir, testPos
  DOUBLE PRECISION, DIMENSION(3)        :: boundPos

  next_cell = 0
  corner = dyn_cell(cell_numb)%corner
  width = dyn_cell(cell_numb)%width
  phot_pos = package(pack_index)%pos
  dir = package(pack_index)%dir
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
  ! we have to choose tTest small enough to be in next cell
  tTest = dist+minwidth/1.D6
  DO
   IF((t1 > 0) .AND. (dist + tTest > t1) .AND. (t1 /= dist)) THEN 
    tTest = tTest/4.E0
   ELSE IF((t2 > 0) .AND. (dist + tTest > t2) .AND. (t2 /= dist)) THEN
    tTest = tTest/4.E0
   ELSE IF((t3 > 0) .AND. (dist + tTest > t3) .AND. (t3 /= dist)) THEN 
    tTest = tTest/4.E0
   ELSE IF((t4 > 0) .AND. (dist + tTest > t4) .AND. (t4 /= dist)) THEN 
    tTest = tTest/4.E0
   ELSE IF((t5 > 0) .AND. (dist + tTest > t5) .AND. (t5 /= dist)) THEN 
    tTest = tTest/4.E0
   ELSE IF((t6 > 0) .AND. (dist + tTest > t6) .AND. (t6 /= dist)) THEN 
    tTest = tTest/4.E0
   ELSE
    EXIT
   END IF
  END DO
  testPos = phot_pos + dir * tTest
  boundPos = phot_pos + dir * dist
  IF((testPos(1) <= -xmax) .OR. (testPos(1) >= xmax)) next_cell = -99
  IF((testPos(2) <= -ymax) .OR. (testPos(2) >= ymax)) next_cell = -99
  IF((testPos(3) <= -zmax) .OR. (testPos(3) >= zmax)) next_cell = -99
  IF(next_cell == 0) CALL find_dyn_cell(testPos,next_cell)
  ! now we will test which part of the cell were hit
  ! walls
!  print*, 'index = ', indx
  END SUBROUTINE next_dyn_cell3
