! this function calculates the spherical coordinates from the cartesian coordinates
! the output is in the form:
!
! r -- the radial distance
! theta -- angle from the z-axe
! phi -- angle from the x-axe (positive direction)
!
! INPUT: cart_coor(DBLE(const_dimofspace)): cartesian coordinates
! OUTPUT: spher_coor(DBLE(const_dimofspace)): spherical coordinates
!
! RETURN: 1x
!
FUNCTION find_spher_coor(cart_coor, spher_coor)

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cart_coor
LOGICAL                                 :: find_spher_coor
DOUBLE PRECISION                        :: coor_x, coor_y, coor_z, coor_r, theta, phi
DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: spher_coor

write(*,*) 'find_spher_coor: start'

coor_x = cart_coor(ind_x)
coor_y = cart_coor(ind_y)
coor_z = cart_coor(ind_z)

coor_r = sqrt(coor_x**2 + coor_y**2 + coor_z**2)
theta = acos(coor_z / sqrt(coor_x**2 + coor_y**2 + coor_z**2))
phi = coor_x / sqrt(coor_x**2 + coor_y**2)

spher_coor(ind_x) = coor_r
spher_coor(ind_y) = theta
spher_coor(ind_z) = phi

find_spher_coor = .TRUE.

! RETURN point
RETURN


END FUNCTION find_spher_coor
