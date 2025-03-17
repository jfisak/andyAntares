SUBROUTINE packet_dynamics(pack_index)


USE types
USE counters

IMPLICIT NONE

INTEGER                                 :: pack_index
INTEGER                                 :: pack_type

pack_type = package(pack_index)%typ

DO  WHILE (package(pack_index)%active == 1)
 IF ((pack_type == type_rpkt) .OR. (pack_type == type_vrpkt)) THEN
  IF(package(pack_index)%n_interactions .EQ. 1000000) THEN
   write(*,*) 'package ', pack_index, ' interacted for 2000000 times and will be destroyed...'
   package(pack_index)%active = 0
   count_des_inte = count_des_inte + 1
  END IF
   ! If the packet is of type rpkt, it represents a photon. So it needs to be propagated.
   CALL do_rpackage(pack_index)
 ELSE IF ((pack_type .EQ. type_kpkt) .OR. (pack_type == type_vkpkt)) THEN 
   ! If the packet is of type kpkt, it represents thermal kinetic energy.
   ! Sample all possible cooling processes and randomly select one of them
    CALL do_kpackage(pack_index)
 ELSE IF ((pack_type .EQ. type_ipkt) .OR. (pack_type == type_vipkt)) THEN 
   ! If the packet is of type ipkt, it represents atomic internal energy (excitation/ionization).
   ! Calculate all transition probabilities and randomly select one of them (macro-atom formalism)
   CALL do_ipackage(pack_index)
 ELSE IF ((pack_type == type_dpkt) .OR. (pack_type == type_vdpkt)) THEN
  ! this state corresponds to diffusive approximation
  CALL do_dpackage(pack_index)
 ELSE
    STOP 'ERROR unknown package typ'
 END IF
END DO


END SUBROUTINE packet_dynamics
