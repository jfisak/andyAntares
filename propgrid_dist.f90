! calculates distance to the boundary of the propGrid from any outer point and direction
!
! INPUT: cur_paket(INT) -- index of packet
! OUTPUT: bound_dist(DBLE) -- distance to the boundary
!
SUBROUTINE propgrid_dist(cur_paket, bound_dist)

USE types
IMPLICIT NONE

TYPE(photon)                                    :: cur_paket
DOUBLE PRECISION                                :: bound_dist
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cross_pos, cur_pos, direction
DOUBLE PRECISION, PARAMETER                     :: largeNumber = 1.D99
DOUBLE PRECISION                                :: txm, txp, typ, tym, tzm, tzp
INTEGER                                         :: next_cross
LOGICAL                                         :: snapped

cur_pos = cur_paket%pos
direction = cur_paket%dir
txp = 0.D0
txm = 0.D0
typ = 0.D0
tym = 0.D0
tzp = 0.D0
tzm = 0.D0

txm = (cur_pos(ind_x) - xmin)/direction(ind_x)
txp = (cur_pos(ind_x) - xmax)/direction(ind_x)
tym = (cur_pos(ind_y) - ymin)/direction(ind_y)
typ = (cur_pos(ind_y) - ymax)/direction(ind_y)
tzm = (cur_pos(ind_z) - zmin)/direction(ind_z)
tzp = (cur_pos(ind_z) - zmax)/direction(ind_z)

if(txm < bound_dist .and. txm > 0.D0) then
 next_cross = negx
 bound_dist = txm
end if
if(txp < bound_dist .and. txp > 0.D0) then
 next_cross = posx
 bound_dist = txp
end if
if(tym < bound_dist .and. tym > 0.D0) then
 next_cross = negy
 bound_dist = tym
end if
if(typ < bound_dist .and. typ > 0.D0) then
 next_cross = posy
 bound_dist = typ
end if
if(tzm < bound_dist .and. txm > 0.D0) then
 next_cross = negz
 bound_dist = tzm
end if
if(tzp < bound_dist .and. txp > 0.D0) then
 next_cross = posz
 bound_dist = tzp
end if

cross_pos = cur_pos + direction * bound_dist
if(next_cross == txm .or. next_cross == txp) then
 ! 
 if(cross_pos(ind_y) > ymin .and. cross_pos(ind_y) < ymax .and. &
  cross_pos(ind_z) > zmin .and. cross_pos(ind_z) < zmax) then
   snapped = .true.
 end if
end if

if(next_cross == tym .or. next_cross == typ) then
 if(cross_pos(ind_x) > xmin .and. cross_pos(ind_x) < xmax .and. &
  cross_pos(ind_z) > zmin .and. cross_pos(ind_z) < zmax) then
   snapped = .true.
 end if
 
end if
if(next_cross == tzm .or. next_cross == tzp) then
 if( cross_pos(ind_y) > ymin .and. cross_pos(ind_y) < ymax .and. &
  cross_pos(ind_x) > xmin .and. cross_pos(ind_x) < xmax) then
   snapped = .true.
 end if
 
end if

if(snapped) then
 package(ind_x)%pos = package(ind_x)%pos + bound_dist * package(ind_x)%dir 
else
end if

END SUBROUTINE propgrid_dist
