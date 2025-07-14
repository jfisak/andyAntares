! calculates the distance to the closest propGrid cell boundary which a packet
! firstly crosses, there are implemented severral numerical corrections in this
! sbr: a correction againts double crossing the same boundary (without scattering)
! and other numerical problems
! except a calculation of the distance it also sets a forbidden boundary of the packet
! in the case of problems with the packet propagation, there are several diagnostics numbers:
! it analyzes distances and counts, how many are positive, negative and zero, based on this
! analysis, the sbr chooses a treatment
! 
! INPUT: pack_index, INT -- index of a packet
!        cell_numb, INT -- index of the current propGrid cell
! OUTPUT: dist, DBLE -- calculated distance
!
SUBROUTINE bound_dist(pack_index, cell_numb, dist, n_pos, n_neg, n_zer, n_par)

USE types
USE constants

IMPLICIT NONE


! input variables
INTEGER                         :: pack_index, cell_numb
! output variable
DOUBLE PRECISION                :: dist
! variables for dynamic cells
DOUBLE PRECISION                :: tnegx, tnegy, tnegz, tposx, tposy, tposz
! perpendicular distances
DOUBLE PRECISION                :: tpnegx, tpnegy, tpnegz, tpposx, tpposy, tpposz
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: corner, upcorner
DOUBLE PRECISION, DIMENSION(const_dimofspace)  :: dir, pos
INTEGER                         :: forbidden
DOUBLE PRECISION, PARAMETER     :: minie = 1e1
DOUBLE PRECISION, PARAMETER     :: epsilon0 = 1e-3

INTEGER                         :: n_pos, n_neg, n_zer, n_par
INTEGER                         :: next_cross, old_cross

DOUBLE PRECISION, PARAMETER     :: velkeCislo = 1.D99
DOUBLE PRECISION                :: mindist
INTEGER                         :: pom_pgi

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
upcorner = dyn_cell(cell_numb)%upcorner
dir = package(pack_index)%dir
pos = package(pack_index)%pos
forbidden = -99
old_cross = package(pack_index)%next_cross
IF(old_cross == posx) THEN
 forbidden = negx
ELSE IF(old_cross == negx) THEN
 forbidden = posx
ELSE IF(old_cross == posy) THEN
 forbidden = negy
ELSE IF(old_cross == negy) THEN
 forbidden = posy
ELSE IF(old_cross == posz) THEN
 forbidden = negz
ELSE IF(old_cross == negz) THEN
 forbidden = posz
END IF
IF(debug == 2) THEN
 CALL find_dyn_cell1(pos, pom_pgi)
 write(*,*) 'bound_dist: cur_pgi = ', pom_pgi
 write(*,*) 'bound_dist: cell starting = ', corner
 write(*,*) 'bound_dist: packet pos = ', pos
 write(*,*) 'bound_dist: cell ending = ', upcorner
END IF


! we will calculate parameters tnegx, ..., tposz
IF(abs(dir(ind_x)) > epsilon0) THEN
 tnegx = (corner(ind_x) - pos(ind_x))/(dir(ind_x))
 tposx = (upcorner(ind_x) - pos(ind_x))/(dir(ind_x))
 IF(tnegx > 0) n_pos = n_pos + 1
 IF(tposx > 0) n_pos = n_pos + 1
 IF(tnegx < 0) n_neg = n_neg + 1
 IF(tposx < 0) n_neg = n_neg + 1
 IF(tnegx == 0.0) n_zer = n_zer + 1
 IF(tposx == 0.0) n_zer = n_zer + 1
ELSE
 package(pack_index)%dir(ind_x) = 0.D0
 tnegx = velkeCislo
 tposx = -velkeCislo
 n_par = n_par + 2
END IF
IF(abs(dir(ind_y)) > epsilon0) THEN
 tnegy = (corner(ind_y) - pos(ind_y))/(dir(ind_y))
 tposy = (upcorner(ind_y) - pos(ind_y))/(dir(ind_y))
 IF(tnegy > 0) n_pos = n_pos + 1
 IF(tposy > 0) n_pos = n_pos + 1
 IF(tnegy < 0) n_neg = n_neg + 1
 IF(tposy < 0) n_neg = n_neg + 1
 IF(tnegy == 0.0) n_zer = n_zer + 1
 IF(tposy == 0.0) n_zer = n_zer + 1
ELSE
 package(pack_index)%dir(ind_y) = 0.D0
 tnegy = velkeCislo
 tposy = -velkeCislo
 n_par = n_par + 2
END IF
IF(abs(dir(ind_z)) > epsilon0) THEN
 tnegz = (corner(ind_z) - pos(ind_z))/(dir(ind_z))
 tposz = (upcorner(ind_z) - pos(ind_z))/(dir(ind_z))
 IF(tnegz > 0) n_pos = n_pos + 1
 IF(tposz > 0) n_pos = n_pos + 1
 IF(tnegz < 0) n_neg = n_neg + 1
 IF(tposz < 0) n_neg = n_neg + 1
 IF(tnegz == 0.0) n_zer = n_zer + 1
 IF(tposz == 0.0) n_zer = n_zer + 1
ELSE
 package(pack_index)%dir(ind_z) = 0.D0
 tnegz = velkeCislo
 tposz = -velkeCislo
 n_par = n_par + 2
END IF

dist = velkeCislo

! we are looking for the bound in front of the photon,
! so we have to choose solution with t > 0
IF( (tnegx > 0.e0) .AND. (tnegx < dist) .AND. forbidden /= negx) THEN
 dist = tnegx
 package(pack_index)%next_cross = negx
END IF
IF( (tnegy > 0.e0)  .AND. (tnegy < dist) .AND. forbidden /= negy) THEN
 dist = tnegy
 package(pack_index)%next_cross = negy
END IF
IF( (tnegz > 0.e0) .AND. (tnegz < dist) .AND. forbidden /=  negz) THEN
 dist = tnegz
 package(pack_index)%next_cross = negz
END IF
IF( (tposx > 0.e0) .AND. (tposx < dist) .AND. forbidden /=  posx) THEN
 dist = tposx
 package(pack_index)%next_cross = posx
END IF
IF( (tposy > 0.e0) .AND. (tposy < dist) .AND. forbidden /=  posy) THEN
 dist = tposy
 package(pack_index)%next_cross = posy
END IF
IF( (tposz > 0.e0) .AND. (tposz < dist) .AND. forbidden /= posz) THEN
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
 tpnegx = pos(ind_x) - corner(ind_x)
 if(abs(tpnegx) < dist) then
  dist = abs(tnegx)
  package(pack_index)%next_cross = negx
 end if
 tpposx = pos(ind_x) - upcorner(ind_x)
 if(abs(tpposx) < dist) then
  dist = abs(tposx)
  package(pack_index)%next_cross = posx
 end if
 tpnegy = pos(ind_y) - corner(ind_y)
 if(abs(tpnegy) < dist) then
  dist = abs(tnegy)
  package(pack_index)%next_cross = negy
 end if
 tpposy = pos(ind_y) - upcorner(ind_y)
 if(abs(tpposy) < dist) then
  dist = abs(tposy)
  package(pack_index)%next_cross = posy
 end if
 tpnegz = pos(ind_z) - corner(ind_z)
 if(abs(tpnegz) < dist) then
  dist = abs(tnegz)
  package(pack_index)%next_cross = negz
 end if
 tpposz = pos(ind_z) - upcorner(ind_z)
 if(abs(tpposz) < dist) then
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

! write(*,*) '*********************************************************************'
! write(*,*) '*********************************************************************'
! write(*,*) 'bound_dist: n_pos = ', n_pos, ' n_neg = ', n_neg
! write(*,*) 'bound_dist: n_zer = ', n_zer, ' n_par = ', n_par
! write(*,*) 'bound_dist: width = ', width
! write(*,*) 'bound_dist: tnegx = ', tnegx, ' tnegy = ', tnegy, ' tnegz = ', tnegz, &
!   ' tposx = ', tposx, ' tposy = ', tposy, 'tposz = ', tposz
! write(*,*) 'bound_dist: tnegx = ', tnegx/width(ind_x), ' tnegy = ', tnegy/width(ind_y), &
! ' tnegz = ', tnegz/width(ind_z), &
!   ' tposx = ', tposx/width(ind_x), ' tposy = ', tposy/width(ind_y), 'tposz = ', tposz/width(ind_z)
! write(*,*) 'bound_dist: dir = ', package(pack_index)%dir
! write(*,*) '*********************************************************************'
! write(*,*) '*********************************************************************'
! write(*,*) 'bound_dist: next_cross = ', package(pack_index)%next_cross
! write(*,*) 'bound_dist: dist = ', dist, ' forbidden = ', forbidden
IF(debug == 2) THEN
 write(*,*) '*********************************************************************'
 write(*,*) '*********************************************************************'
 write(*,*) 'bound_dist: n_pos = ', n_pos, ' n_neg = ', n_neg
 write(*,*) 'bound_dist: n_zer = ', n_zer, ' n_par = ', n_par
 write(*,*) 'bound_dist: tnegx = ', tnegx, ' tnegy = ', tnegy, ' tnegz = ', tnegz, &
   ' tposx = ', tposx, ' tposy = ', tposy, 'tposz = ', tposz
 write(*,*) 'bound_dist: dir = ', package(pack_index)%dir
 write(*,*) '*********************************************************************'
 write(*,*) '*********************************************************************'
 write(*,*) 'bound_dist: next_cross = ', package(pack_index)%next_cross
 write(*,*) 'bound_dist: dist = ', dist, ' forbidden = ', forbidden
END IF


END SUBROUTINE bound_dist
