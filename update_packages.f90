! goes throught all packets and calls sbr packet_dynamics which propagates them through the computational
! domain, also saves temporary packets
!
! INPUT: n_pack(INT) -- number of packets
! OUTPUT: NONE
!
SUBROUTINE update_packages(n_pack)

USE types

IMPLICIT NONE    

INTEGER             :: n_pack, pack_index
  
DO pack_index = tot_saved_packets + 1, n_pack
 ! IF (MODULO(pack_index,100000) .EQ. 0) write(99,*) 'Working on packet ', pack_index,' ...'
 IF (MODULO(pack_index,100000) .EQ. 0) write(*,*) 'Working on packet ', pack_index,' ...'
 ! write(*,*) 'Working on packet ', pack_index,' ...'
 CALL packet_dynamics(pack_index)
 ! Do this loop until something happened with package
 ! saving temporary packet data
 IF(MODULO(pack_index, n_pack_save) == 0) THEN
  CALL save_temp_packs(pack_index)
 END IF
END DO
  
END SUBROUTINE update_packages
