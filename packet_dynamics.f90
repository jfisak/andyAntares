SUBROUTINE packet_dynamics(pack_index)


USE types
USE counters

IMPLICIT NONE

INTEGER                                 :: pack_index
INTEGER                                 :: pack_type
INTEGER, PARAMETER                      :: min_n_interact = 2000000
INTEGER                                 :: testing_counter
INTEGER, PARAMETER                      :: n_testingloop = 2000000
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_pos, old_pos
DOUBLE PRECISION                                :: cur_dist
DOUBLE PRECISION, PARAMETER                     :: min_dist = 1.D1
LOGICAL                                         :: procout=.false.

pack_type = package(pack_index)%typ

testing_counter = 0

DO  WHILE (package(pack_index)%active == 1)
 pack_type = package(pack_index)%typ
 testing_counter = testing_counter + 1
 IF(testing_counter == 1) THEN
  old_pos = package(pack_index)%pos
 END IF
 IF((pack_type == type_rpkt) .OR. (pack_type == type_vrpkt)) THEN
  IF(package(pack_index)%n_interactions .EQ. min_n_interact) THEN
   write(*,*) 'package ', pack_index, ' interacted for ', min_n_interact, ' times and will be destroyed...'
   package(pack_index)%active = 0
   count_des_inte = count_des_inte + 1
  END IF
  ! If the packet is of type rpkt, it represents a photon. So it needs to be propagated.
  if(procout) write(*,*) 'packet_dynamics: r-packet'
  CALL do_rpackage(pack_index)
 ELSE IF ((pack_type .EQ. type_kpkt) .OR. (pack_type == type_vkpkt)) THEN 
  ! If the packet is of type kpkt, it represents thermal kinetic energy.
  ! Sample all possible cooling processes and randomly select one of them
  if(procout) write(*,*) 'packet_dynamics: k-packet'
  CALL do_kpackage(pack_index)
 ELSE IF ((pack_type .EQ. type_ipkt) .OR. (pack_type == type_vipkt)) THEN 
  ! If the packet is of type ipkt, it represents atomic internal energy (excitation/ionization).
  ! Calculate all transition probabilities and randomly select one of them (macro-atom formalism)
  if(procout) write(*,*) 'packet_dynamics: i-packet'
  CALL do_ipackage(pack_index)
 ELSE IF ((pack_type == type_dpkt) .OR. (pack_type == type_vdpkt)) THEN
  ! this state corresponds to diffusive approximation
  if(procout) write(*,*) 'packet_dynamics: d-packet'
  ! write(*,*) 'packet_dynamics: calling do_dpackage for pack_index = ', pack_index
  CALL do_dpackage(pack_index)
 ELSE
    STOP 'ERROR unknown package typ'
 END IF
 IF(testing_counter == n_testingloop) THEN
  cur_pos = package(pack_index)%pos
  cur_dist = sqrt((cur_pos(ind_x) - old_pos(ind_x))**2 + (cur_pos(ind_y) - old_pos(ind_y))**2 + &
   (cur_pos(ind_z) - old_pos(ind_z))**2)
  IF(cur_dist < min_dist .and. pack_type /= type_dpkt) THEN
   write(*,*) 'packet_dynamics: packet ', pack_index, ' is stucked'
   STOP
  END IF
 END IF
END DO


END SUBROUTINE packet_dynamics
