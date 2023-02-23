SUBROUTINE bound_dist(pack_index, cell_numb, dist)

USE types
USE constants

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
DOUBLE PRECISION, PARAMETER     :: minie = 1e1

INTEGER                         :: n_pos, n_neg, n_zer, n_par

DOUBLE PRECISION, PARAMETER     :: velkeCislo = 1.D99
DOUBLE PRECISION                :: mindist

t1 = 0.E0
t2 = 0.E0
t3 = 0.E0
t4 = 0.E0
t5 = 0.E0
t6 = 0.E0

n_pos = 0
n_neg = 0
n_zer = 0
n_par = 0

corner = dyn_cell(cell_numb)%corner
width = dyn_cell(cell_numb)%width
dir = package(pack_index)%dir
pos = package(pack_index)%pos
forbidden = package(pack_index)%next_cross

! we will calculate parameters t1, ..., t6
IF(dir(1) /= 0) THEN
 t1 = (corner(1) - pos(1))/(dir(1))
 t4 = (corner(1) + width(1) - pos(1))/(dir(1))
 IF(t1 > 0) n_pos = n_pos + 1
 IF(t4 > 0) n_pos = n_pos + 1
 IF(t1 < 0) n_neg = n_neg + 1
 IF(t4 < 0) n_neg = n_neg + 1
 IF(t1 == 0.0) n_zer = n_zer + 1
 IF(t4 == 0.0) n_zer = n_zer + 1
ELSE
 t1 = velkeCislo
 t4 = -velkeCislo
 n_par = n_par + 2
END IF
IF(dir(2) /= 0) THEN
 t2 = (corner(2) - pos(2))/(dir(2))
 t5 = (corner(2) + width(2) - pos(2))/(dir(2))
 IF(t2 > 0) n_pos = n_pos + 1
 IF(t5 > 0) n_pos = n_pos + 1
 IF(t2 < 0) n_neg = n_neg + 1
 IF(t5 < 0) n_neg = n_neg + 1
 IF(t2 == 0.0) n_zer = n_zer + 1
 IF(t5 == 0.0) n_zer = n_zer + 1
ELSE
 t2 = velkeCislo
 t5 = -velkeCislo
 n_par = n_par + 2
END IF
IF(dir(3) /= 0) THEN
 t3 = (corner(3) - pos(3))/(dir(3))
 t6 = (corner(3) + width(3) - pos(3))/(dir(3))
 IF(t3 > 0) n_pos = n_pos + 1
 IF(t6 > 0) n_pos = n_pos + 1
 IF(t3 < 0) n_neg = n_neg + 1
 IF(t6 < 0) n_neg = n_neg + 1
 IF(t3 == 0.0) n_zer = n_zer + 1
 IF(t6 == 0.0) n_zer = n_zer + 1
ELSE
 t3 = velkeCislo
 t6 = -velkeCislo
 n_par = n_par + 2
END IF


dist = velkeCislo

! we are looking for the bound in front of the photon,
! so we have to choose solution with t > 0
IF( (t1 > 0.e0) .AND. (t1 < dist) .AND. forbidden /= posx) THEN
 dist = t1
 package(pack_index)%next_cross = negx
END IF
IF( (t2 > 0.e0)  .AND. (t2 < dist)  .AND. forbidden /= posy) THEN
 dist = t2
 package(pack_index)%next_cross = negy
END IF
IF( (t3 > 0.e0) .AND. (t3 < dist)  .AND. forbidden /=  posz) THEN
 dist = t3
 package(pack_index)%next_cross = negz
END IF
IF( (t4 > 0.e0) .AND. (t4 < dist)  .AND. forbidden /=  negx) THEN
 dist = t4
 package(pack_index)%next_cross = posx
END IF
IF( (t5 > 0.e0) .AND. (t5 < dist)  .AND. forbidden /=  negy) THEN
 dist = t5
 package(pack_index)%next_cross = posy
END IF
IF( (t6 > 0.e0) .AND. (t6 < dist)  .AND. forbidden /= negz) THEN
 dist = t6
 package(pack_index)%next_cross = posz
END IF

IF(forbidden == -99 .and. n_pos > 3) THEN
 mindist = dist
 dist = velkeCislo
 IF( (t1 > mindist) .AND. (t1 < dist)) THEN
  dist = t1
  package(pack_index)%next_cross = negx
 END IF
 IF( (t2 > mindist)  .AND. (t2 < dist)) THEN
  dist = t2
  package(pack_index)%next_cross = negy
 END IF
 IF( (t3 > mindist) .AND. (t3 < dist)) THEN
  dist = t3
  package(pack_index)%next_cross = negz
 END IF
 IF( (t4 > mindist) .AND. (t4 < dist)) THEN
  dist = t4
  package(pack_index)%next_cross = posx
 END IF
 IF( (t5 > mindist) .AND. (t5 < dist)) THEN
  dist = t5
  package(pack_index)%next_cross = posy
 END IF
 IF( (t6 > mindist) .AND. (t6 < dist)) THEN
  dist = t6
  package(pack_index)%next_cross = posz
 END IF
END IF

IF(n_neg > 3) THEN
 dist = velkeCislo
 ! calculation of perpendicular distance to propGrid cell surfaces
 t1 = pos(1) - corner(1)
 if(abs(t1) < dist) then
  dist = abs(t1)
  package(pack_index)%next_cross = negx
 end if
 t4 = pos(1) - corner(1) - width(1)
 if(abs(t4) < dist) then
  dist = abs(t4)
  package(pack_index)%next_cross = posx
 end if
 t2 = pos(2) - corner(2)
 if(abs(t2) < dist) then
  dist = abs(t2)
  package(pack_index)%next_cross = negy
 end if
 t5 = pos(2) - corner(2) - width(2)
 if(abs(t5) < dist) then
  dist = abs(t5)
  package(pack_index)%next_cross = posy
 end if
 t3 = pos(3) - corner(3)
 if(abs(t3) < dist) then
  dist = abs(t3)
  package(pack_index)%next_cross = negz
 end if
 t6 = pos(3) - corner(3) - width(3)
 if(abs(t6) < dist) then
  dist = abs(t6)
  package(pack_index)%next_cross = posz
 end if
 dist = -dist
END IF

IF(n_zer == 1 .and. n_pos < 3 .and. n_zer == 0) THEN
 dist = -1.0
 IF(t1 == 0.0) package(pack_index)%next_cross = negx
 IF(t2 == 0.0) package(pack_index)%next_cross = negy
 IF(t3 == 0.0) package(pack_index)%next_cross = negz
 IF(t4 == 0.0) package(pack_index)%next_cross = posx
 IF(t5 == 0.0) package(pack_index)%next_cross = posy
 IF(t6 == 0.0) package(pack_index)%next_cross = posz
END IF

IF(debug == 2) THEN
 write(*,*) 'bound_dist: n_pos = ', n_pos, ' n_neg = ', n_neg
 write(*,*) 'bound_dist: t1 = ', t1, ' t2 = ', t2, ' t3 = ', t3, ' t4 = ', t4, ' t5 = ', t5, ' t6 = ', t6
 write(*,*) 'bound_dist: dist = ', dist
END IF


END SUBROUTINE bound_dist
