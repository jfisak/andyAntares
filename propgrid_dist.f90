SUBROUTINE propgrid_dist(cur_paket, bound_dist)

USE types
IMPLICIT NONE

TYPE(photon)                                    :: cur_paket
DOUBLE PRECISION                                :: bound_dist
DOUBLE PRECISION, DIMENSION(3)                  :: cross_pos, cur_pos, direction
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

txm = (cur_pos(1) - xmin)/direction(1)
txp = (cur_pos(1) - xmax)/direction(1)
tym = (cur_pos(2) - ymin)/direction(2)
typ = (cur_pos(2) - ymax)/direction(2)
tzm = (cur_pos(3) - zmin)/direction(3)
tzp = (cur_pos(3) - zmax)/direction(3)

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
 if(cross_pos(2) > ymin .and. cross_pos(2) < ymax .and. &
  cross_pos(3) > zmin .and. cross_pos(3) < zmax) then
   snapped = .true.
 end if
end if

if(next_cross == tym .or. next_cross == typ) then
 if(cross_pos(1) > xmin .and. cross_pos(1) < xmax .and. &
  cross_pos(3) > zmin .and. cross_pos(3) < zmax) then
   snapped = .true.
 end if
 
end if
if(next_cross == tzm .or. next_cross == tzp) then
 if( cross_pos(2) > ymin .and. cross_pos(2) < ymax .and. &
  cross_pos(1) > xmin .and. cross_pos(1) < xmax) then
   snapped = .true.
 end if
 
end if

if(snapped) then
 package(1)%pos = package(1)%pos + bound_dist * package(1)%dir 
else
end if

END SUBROUTINE propgrid_dist
