SUBROUTINE find_dist(pack_index, cell_numb, dist)

USE types

IMPLICIT NONE


    ! input variables
    INTEGER                         :: pack_index, cell_numb
    ! output variable
    DOUBLE PRECISION                :: dist
    ! variables for dynamic cells
    DOUBLE PRECISION                :: t1, t2, t3, t4, t5, t6
    DOUBLE PRECISION, DIMENSION(3)  :: corner, width
    DOUBLE PRECISION, DIMENSION(3)  :: dir

  corner = dyn_cell(cell_numb)%corner
  width = dyn_cell(cell_numb)%width
  cell_numb = package(pack_index)%cell_numb
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

END SUBROUTINE find_dist
