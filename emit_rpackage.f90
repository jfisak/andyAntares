SUBROUTINE emit_rpackage(pack_index)

USE types

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: D
DOUBLE PRECISION, DIMENSION(3)    :: cmf_direction, rf_direction, vel_vec

DOUBLE PRECISION, DIMENSION(3)          :: cur_pos
INTEGER                                 :: cur_mgi, get_package_model_index
  
cur_pos = package(pack_index)%pos
cur_mgi = get_package_model_index(pack_index)

! Make it an r-pkt and make sure that all cell transitions are
! possible therafter
package(pack_index)%typ = type_rpkt
package(pack_index)%next_cross = NONE

! Randomly emitt in CMF
CALL random_unitvector(cmf_direction) ! Only for cmf

! Transform direction from cmf to rf. Take negativ velocity for the
! aberration formula (rf to cmf trafo would require positive
! velocity). See e.g. Mihalas and Mihalas Eq. 89.6
CALL velo(pack_index, cur_pos, cur_mgi, vel_vec, velApprox)
CALL angle_aberration(cmf_direction, -1.D0*vel_vec, rf_direction)
package(pack_index)%dir = rf_direction

! Finally update the packets frequency and energy
! See e.g. Mihalas and Mihalas Eq. 89.5
CALL doppler_factor(pack_index, D)

package(pack_index)%e_rf = package(pack_index)%e_cmf / D
package(pack_index)%freq_rf =  package(pack_index)%freq_cmf / D

END SUBROUTINE emit_rpackage
