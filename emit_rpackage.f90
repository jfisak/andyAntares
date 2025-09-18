! emit a new r-packet in a random direction
!
! INPUT: pack_index(INT): index of a packet
! OUTPUT: NONE
!
SUBROUTINE emit_rpackage(pack_index)

USE types
USE constants
USE dummypacket

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: doppler_D
DOUBLE PRECISION, DIMENSION(const_dimofspace)    :: cmf_direction, rf_direction, vel_vec

DOUBLE PRECISION, DIMENSION(const_dimofspace)          :: cur_pos
INTEGER                                 :: dummy_pack_index
INTEGER                                 :: cur_mgi, get_package_model_index
  
IF(pack_index <= SIZE(package)) THEN
 cur_pos = package(pack_index)%pos
 cur_mgi = get_package_model_index(pack_index)
 ! Make it an r-pkt and make sure that all cell transitions are
 ! possible therafter
 package(pack_index)%typ = type_rpkt
 package(pack_index)%next_cross = NONE
ELSE
 dummy_pack_index = pack_index - SIZE(package)
 cur_pos = dummypackage(dummy_pack_index)%pos
 cur_mgi = get_package_model_index(dummy_pack_index)
 ! Make it an r-pkt and make sure that all cell transitions are
 ! possible therafter
 dummypackage(dummy_pack_index)%typ = type_rpkt
 dummypackage(dummy_pack_index)%next_cross = NONE
END IF


! Randomly emitt in CMF
CALL random_unitvector(cmf_direction) ! Only for cmf

! Transform direction from cmf to rf. Take negativ velocity for the
! aberration formula (rf to cmf trafo would require positive
! velocity). See e.g. Mihalas and Mihalas Eq. 89.6
CALL velo(pack_index, vel_vec, velApprox)
CALL angle_aberration(cmf_direction, -1.D0*vel_vec, rf_direction)
IF(pack_index <= SIZE(package)) THEN
 package(pack_index)%dir = rf_direction
ELSE
 dummypackage(dummy_pack_index)%dir = rf_direction
END IF

! Finally update the packets frequency and energy
! See e.g. Mihalas and Mihalas Eq. 89.5
CALL doppler_factor(pack_index, doppler_D)

IF(pack_index <= SIZE(package)) THEN
 package(pack_index)%e_rf = package(pack_index)%e_cmf / doppler_D
 package(pack_index)%freq_rf =  package(pack_index)%freq_cmf / doppler_D
ELSE
 dummypackage(dummy_pack_index)%e_rf = package(dummy_pack_index)%e_cmf / doppler_D
 dummypackage(dummy_pack_index)%freq_rf =  package(dummy_pack_index)%freq_cmf / doppler_D
END IF

END SUBROUTINE emit_rpackage
