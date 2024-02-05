SUBROUTINE brtm()

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: obs_point, ccd_point
INTEGER                                                 :: cur_vpack
INTEGER, PARAMETER                                      :: Nvpackets = 100

DEALLOCATE(package)

! observing point
obs_point = (/ -2*R_inf  ,  0.D0,  0.D0 /)
                                    
ccd_point = (/ -1.9*R_inf,  0.D0,  0.D0 /)


ALLOCATE(package(Nvpackets))
! CCD chip


cur_vpack = 1

! a basic initialisation of a packet
package(cur_vpack)%pos = obs_point
package(cur_vpack)%dir = (ccd_point - obs_point)/norm2(ccd_point - obs_point)

! basic properties of a packet
package(cur_vpack)%active = 1
package(cur_vpack)%typ = type_rpkt
package(cur_vpack)%n_interactions = 0
package(cur_vpack)%next_cross = NONE
package(cur_vpack)%e_rf = 1.e10
package(cur_vpack)%last_line = no_line
package(cur_vpack)%delta_s = 0.D0

! calling a sbr to process a v-packet
CALL do_vpackage(cur_vpack)

END SUBROUTINE brtm
