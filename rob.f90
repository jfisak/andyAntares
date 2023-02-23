! red or blue shift
FUNCTION rob(pack_index, b_dist)

USE types
USE constants
IMPLICIT NONE


INTEGER                         :: pack_index
DOUBLE PRECISION                :: b_dist

INTEGER                         :: dummypackage

LOGICAL                         :: redshift, rob

DOUBLE PRECISION                :: nu_0, nu_b

dummypackage = SIZE(package)

package(dummypackage) = package(pack_index)

nu_0 = package(pack_index)%freq_cmf

CALL move_package(dummypackage, b_dist, 0, .false.)

nu_b = package(dummypackage)%freq_cmf

! write(*,*) 'rob: nu_0 = ', nu_0, ' nu_b = ', nu_b

if(nu_b < nu_0) THEN
 redshift = .TRUE.
else
 redshift = .FALSE.
end if

rob = redshift
RETURN

END FUNCTION rob
