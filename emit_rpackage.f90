  SUBROUTINE emit_rpackage(pack_index)

  USE types

  IMPLICIT NONE    

    INTEGER                           :: pack_index
    DOUBLE PRECISION                  :: D
    DOUBLE PRECISION, DIMENSION(3)    :: cmf_direction, rf_direction, vel_vec

!Make it an r-pkt and make sure that all cell transitions are possible therafter
    package(pack_index)%typ = type_rpkt
    package(pack_index)%last_cross = NONE


!Randomly emitt in CMF
    CALL random_unitvector(cmf_direction) ! Only for cmf

!   Transform direction from cmf to rf. Take negativ velocity for the aberration formula
!   (rf to cmf trafo would require positive velocity)
    CALL velo(pack_index,vel_vec)
    CALL angle_aberration(cmf_direction, -1.D0*vel_vec, rf_direction)
    package(pack_index)%dir = rf_direction

!   Finally update the packets frequency and energy
    CALL doppler_factor(pack_index, D)
    
    package(pack_index)%e_rf = package(pack_index)%e_cmf / D
    package(pack_index)%freq_rf =  package(pack_index)%freq_cmf / D

    
  END SUBROUTINE emit_rpackage

