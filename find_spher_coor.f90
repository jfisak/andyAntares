FUNCTION find_spher_coor(cart_coor, spher_coor)
! this function calculates the spherical coordinates from the cartesian coordinates
! the output is in the form:
!
! r -- the radial distance
! theta -- angle from the z-axe
! phi -- angle from the x-axe (positive direction)

DOUBLE PRECISION, DIMENSION(3)          :: cart_coor
LOGICAL                                 :: find_spher_coor
DOUBLE PRECISION                        :: x, y, z, r, theta, phi
DOUBLE PRECISION, DIMENSION(3)          :: spher_coor

write(*,*) 'find_spher_coor: start'

x = cart_coor(1)
y = cart_coor(2)
z = cart_coor(3)

r = sqrt(x**2 + y**2 + z**2)
theta = acos(z / sqrt(x**2 + y**2 + z**2))
phi = x / sqrt(x**2 + y**2)

spher_coor(1) = r
spher_coor(2) = theta
spher_coor(3) = phi

find_spher_coor = .TRUE.

RETURN


END FUNCTION find_spher_coor
