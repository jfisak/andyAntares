SUBROUTINE bound_dist(pack_index, cell_numb, dist)

USE types
USE constants

IMPLICIT NONE


! input variables
INTEGER                         :: pack_index, cell_numb
! output variable
DOUBLE PRECISION                :: dist
! variables for dynamic cells
DOUBLE PRECISION                :: tnegx, tnegy, tnegz, tposx, tposy, tposz
DOUBLE PRECISION, DIMENSION(3)  :: corner, width
DOUBLE PRECISION, DIMENSION(3)  :: dir, pos
INTEGER                         :: forbidden
DOUBLE PRECISION, PARAMETER     :: minie = 1e1

INTEGER                         :: n_pos, n_neg, n_zer, n_par
INTEGER                         :: next_cross

DOUBLE PRECISION, PARAMETER     :: velkeCislo = 1.D99
DOUBLE PRECISION                :: mindist

tnegx = 0.E0
tnegy = 0.E0
tnegz = 0.E0
tposx = 0.E0
tposy = 0.E0
tposz = 0.E0

n_pos = 0
n_neg = 0
n_zer = 0
n_par = 0

corner = dyn_cell(cell_numb)%corner
width = dyn_cell(cell_numb)%width
dir = package(pack_index)%dir
pos = package(pack_index)%pos
forbidden = package(pack_index)%next_cross

! we will calculate parameters tnegx, ..., tposz
IF(dir(1) /= 0) THEN
 tnegx = (corner(1) - pos(1))/(dir(1))
 tposx = (corner(1) + width(1) - pos(1))/(dir(1))
 IF(tnegx > 0) n_pos = n_pos + 1
 IF(tposx > 0) n_pos = n_pos + 1
 IF(tnegx < 0) n_neg = n_neg + 1
 IF(tposx < 0) n_neg = n_neg + 1
 IF(tnegx == 0.0) n_zer = n_zer + 1
 IF(tposx == 0.0) n_zer = n_zer + 1
ELSE
 tnegx = velkeCislo
 tposx = -velkeCislo
 n_par = n_par + 2
END IF
IF(dir(2) /= 0) THEN
 tnegy = (corner(2) - pos(2))/(dir(2))
 tposy = (corner(2) + width(2) - pos(2))/(dir(2))
 IF(tnegy > 0) n_pos = n_pos + 1
 IF(tposy > 0) n_pos = n_pos + 1
 IF(tnegy < 0) n_neg = n_neg + 1
 IF(tposy < 0) n_neg = n_neg + 1
 IF(tnegy == 0.0) n_zer = n_zer + 1
 IF(tposy == 0.0) n_zer = n_zer + 1
ELSE
 tnegy = velkeCislo
 tposy = -velkeCislo
 n_par = n_par + 2
END IF
IF(dir(3) /= 0) THEN
 tnegz = (corner(3) - pos(3))/(dir(3))
 tposz = (corner(3) + width(3) - pos(3))/(dir(3))
 IF(tnegz > 0) n_pos = n_pos + 1
 IF(tposz > 0) n_pos = n_pos + 1
 IF(tnegz < 0) n_neg = n_neg + 1
 IF(tposz < 0) n_neg = n_neg + 1
 IF(tnegz == 0.0) n_zer = n_zer + 1
 IF(tposz == 0.0) n_zer = n_zer + 1
ELSE
 tnegz = velkeCislo
 tposz = -velkeCislo
 n_par = n_par + 2
END IF


dist = velkeCislo

! we are looking for the bound in front of the photon,
! so we have to choose solution with t > 0
IF( (tnegx > 0.e0) .AND. (tnegx < dist) .AND. forbidden /= posx) THEN
 dist = tnegx
 package(pack_index)%next_cross = negx
END IF
IF( (tnegy > 0.e0)  .AND. (tnegy < dist) .AND. forbidden /= posy) THEN
 dist = tnegy
 package(pack_index)%next_cross = negy
END IF
IF( (tnegz > 0.e0) .AND. (tnegz < dist) .AND. forbidden /=  posz) THEN
 dist = tnegz
 package(pack_index)%next_cross = negz
END IF
IF( (tposx > 0.e0) .AND. (tposx < dist) .AND. forbidden /=  negx) THEN
 dist = tposx
 package(pack_index)%next_cross = posx
END IF
IF( (tposy > 0.e0) .AND. (tposy < dist) .AND. forbidden /=  negy) THEN
 dist = tposy
 package(pack_index)%next_cross = posy
END IF
IF( (tposz > 0.e0) .AND. (tposz < dist) .AND. forbidden /= negz) THEN
 dist = tposz
 package(pack_index)%next_cross = posz
END IF

IF(forbidden == -99 .and. n_pos > 3) THEN
 mindist = dist
 dist = velkeCislo
 IF( (tnegx > mindist) .AND. (tnegx < dist)) THEN
  dist = tnegx
  package(pack_index)%next_cross = negx
 END IF
 IF( (tnegy > mindist)  .AND. (tnegy < dist)) THEN
  dist = tnegy
  package(pack_index)%next_cross = negy
 END IF
 IF( (tnegz > mindist) .AND. (tnegz < dist)) THEN
  dist = tnegz
  package(pack_index)%next_cross = negz
 END IF
 IF( (tposx > mindist) .AND. (tposx < dist)) THEN
  dist = tposx
  package(pack_index)%next_cross = posx
 END IF
 IF( (tposy > mindist) .AND. (tposy < dist)) THEN
  dist = tposy
  package(pack_index)%next_cross = posy
 END IF
 IF( (tposz > mindist) .AND. (tposz < dist)) THEN
  dist = tposz
  package(pack_index)%next_cross = posz
 END IF
END IF


! this part returns negative number if more than three distances are negative
! it means that the packet is located in the neighboring cell and its direction has been
! changed
IF(n_neg > 3) THEN
 dist = velkeCislo
 ! calculation of perpendicular distance to propGrid cell surfaces
 tnegx = pos(1) - corner(1)
 if(abs(tnegx) < dist) then
  dist = abs(tnegx)
  package(pack_index)%next_cross = negx
 end if
 tposx = pos(1) - corner(1) - width(1)
 if(abs(tposx) < dist) then
  dist = abs(tposx)
  package(pack_index)%next_cross = posx
 end if
 tnegy = pos(2) - corner(2)
 if(abs(tnegy) < dist) then
  dist = abs(tnegy)
  package(pack_index)%next_cross = negy
 end if
 tposy = pos(2) - corner(2) - width(2)
 if(abs(tposy) < dist) then
  dist = abs(tposy)
  package(pack_index)%next_cross = posy
 end if
 tnegz = pos(3) - corner(3)
 if(abs(tnegz) < dist) then
  dist = abs(tnegz)
  package(pack_index)%next_cross = negz
 end if
 tposz = pos(3) - corner(3) - width(3)
 if(abs(tposz) < dist) then
  dist = abs(tposz)
  package(pack_index)%next_cross = posz
 end if
 dist = -dist
END IF

! the packet is in the edge or in the corner of the propGrid cells
IF(n_zer > 1) THEN
 IF((tnegy == 0.D0 .or. tposy == 0) .and. (tnegz == 0.D0 .or. tposz == 0)) THEN
  next_cross = edyz
 ELSE IF((tnegx == 0.D0 .or. tposx == 0) .and. (tnegz == 0.D0 .or. tposz == 0)) THEN
  next_cross = edxz
 ELSE IF((tnegx == 0.D0 .or. tposx == 0) .and. (tnegy == 0.D0 .or. tposy == 0)) THEN
  next_cross = edxy
 END IF
 package(pack_index)%next_cross = next_cross
 dist = -1.D0
END IF

IF(debug == 2) THEN
 write(*,*) '*********************************************************************'
 write(*,*) '*********************************************************************'
 write(*,*) 'bound_dist: n_pos = ', n_pos, ' n_neg = ', n_neg
 write(*,*) 'bound_dist: n_zer = ', n_zer, ' n_par = ', n_par
 write(*,*) 'bound_dist: tnegx = ', tnegx/R_star, ' tnegy = ', tnegy/R_star, ' tnegz = ', tnegz/R_star, &
   ' tposx = ', tposx/R_star, ' tposy = ', tposy/R_star, 'tposz = ', tposz/R_star
 write(*,*) '*********************************************************************'
 write(*,*) '*********************************************************************'
 write(*,*) 'bound_dist: next_cross = ', package(pack_index)%next_cross
 write(*,*) 'bound_dist: dist = ', dist/R_star, ' forbidden = ', forbidden
END IF


END SUBROUTINE bound_dist
