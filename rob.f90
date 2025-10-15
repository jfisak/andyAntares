! red or blue shift
! 
! INPUT: pack_index(INT): index of packet
!        b_dist(DBL): bound distance
! OUTPUT: rob(LOG): red shift (T), blue shift (F)
!
! 1x RETURN point
!
FUNCTION rob(pack_index, b_dist)

USE types
USE constants
USE dummypacket
IMPLICIT NONE


INTEGER                         :: pack_index
DOUBLE PRECISION                :: b_dist

INTEGER                         :: dummypack_index, cur_dummypack

LOGICAL                         :: redshift, rob

DOUBLE PRECISION                :: nu_0, nu_b

cur_dummypack = find_free_index()
CALL copy_package(pack_index, cur_dummypack)
dummypack_index = cur_dummypack + SIZE(package)

nu_0 = package(pack_index)%freq_cmf

CALL move_package(dummypack_index, b_dist, 0, .false.)

nu_b = dummypackage(cur_dummypack)%freq_cmf

! write(*,*) 'rob: nu_0 = ', nu_0, ' nu_b = ', nu_b

if(nu_b < nu_0) THEN
 redshift = .TRUE.
else
 redshift = .FALSE.
end if

CALL deactivate_dummy_packet(cur_dummypack)

rob = redshift

! RETURN point
RETURN

END FUNCTION rob
