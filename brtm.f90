SUBROUTINE brtm()

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: obs_point, ccd_point
INTEGER                                                 :: cur_vpack
INTEGER, PARAMETER                                      :: Nvpackets = 100


! observing point
obs_point = (/ 0.D0, -2*R_inf, 0.D0 /)

ccd_point = (/ 0.D0, -1.9*R_inf, 0.D0 /)


ALLOCATE(vpackage(Nvpackets))
! CCD chip




cur_vpack = 1

vpackage(cur_vpack)%active = 1
vpackage(cur_vpack)%pos = obs_point
vpackage(cur_vpack)%dir = (ccd_point - obs_point)/norm2(ccd_point - obs_point)
CALL do_vpackage(cur_vpack)


STOP 'brtm testing'



END SUBROUTINE brtm

