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
    DOUBLE PRECISION, DIMENSION(3)  :: dir, pos
    INTEGER                         :: forbidden

  cell_numb = package(pack_index)%cell_numb
  corner = dyn_cell(cell_numb)%corner
  width = dyn_cell(cell_numb)%width
  dir = package(pack_index)%dir
  pos = package(pack_index)%pos
  forbidden = package(pack_index)%next_cross
 ! we will calculate parameters t1,...,t6
 IF(dir(1) /= 0) THEN
  t1 = (corner(1) - pos(1))/(dir(1))
  t4 = (corner(1) + width(1) - pos(1))/(dir(1))
 ELSE
  t1 = 0
  t4 = 0
 END IF
 IF(dir(2) /= 0) THEN
  t2 = (corner(2) - pos(2))/(dir(2))
  t5 = (corner(2) + width(2) - pos(2))/(dir(2))
 ELSE
  t2 = 0
  t5 = 0
 END IF
 IF(dir(3) /= 0) THEN
  t3 = (corner(3) - pos(3))/(dir(3))
  t6 = (corner(3) + width(3) - pos(3))/(dir(3))
 ELSE
  t3 = 0
  t6 = 0
 END IF

  dist = 1.D99
  ! we are looking for the bound in front of the photon,
  ! so we have to choose solution with t > 0
  IF( (t1 > 0.E0) .AND. (t1 < dist) .AND. forbidden /= posx) THEN
   dist = t1
   package(pack_index)%next_cross = negx
  END IF
  IF( (t2 > 0.E0)  .AND. (t2 < dist)  .AND. forbidden /= posy) THEN
   dist = t2
   package(pack_index)%next_cross = negy
  END IF
  IF( (t3 > 0.E0) .AND. (t3 < dist)  .AND. forbidden /=  posz) THEN
   dist = t3
   package(pack_index)%next_cross = negz
  END IF
  IF( (t4 > 0.E0) .AND. (t4 < dist)  .AND. forbidden /=  negx) THEN
   dist = t4
   package(pack_index)%next_cross = posx
  END IF
  IF( (t5 > 0.E0) .AND. (t5 < dist)  .AND. forbidden /=  negy) THEN
   dist = t5
   package(pack_index)%next_cross = posy
  END IF
  IF( (t6 > 0.E0) .AND. (t6 < dist)  .AND. forbidden /= negz) THEN
   dist = t6
   package(pack_index)%next_cross = posz
  END IF
   !print*, 'find_dist: t1 = ', t1, ' t2 = ', t2, ' t3 = ', t3, ' t4 = ', t4, ' t5 = ', t5, ' t6 = ', t6

END SUBROUTINE find_dist
