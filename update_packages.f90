 SUBROUTINE update_packages(n_pack)

  USE types
USE constants
  USE counters
!  USE rand2

  IMPLICIT NONE    

  INTEGER             :: n_pack, pack_index
  
! OPEN(UNIT=36, FILE='macroatom.dat')
!  OPEN (UNIT=3, FILE='position.dat')
! write(*,*) 'WARNING: positions of packets are being written into the file, if the number of &
!                packages is large the file will be very large'
DO pack_index = tot_saved_packets + 1, n_pack
 ! IF(initrs(my_rank + 1) .EQV. .FALSE.) THEN
 !  CALL init_random_seed()
 !  initrs = .TRUE.
 ! END IF
  ! write(36, *) 'r-packet: ', pack_index
  ! write(*,*) 'update_packages: package = ', pack_index
 IF (MODULO(pack_index,100000) .EQ. 0) write(99,*) 'Working on packet ', pack_index,' ...'
 ! IF (MODULO(pack_index,10) .EQ. 0) write(*,*) 'Working on packet ', pack_index,' ...'
 ! write(*,*) 'Working on packet ', pack_index,' ...'
 CALL packet_dynamics(pack_index)
 ! Do this loop until something happened with package
 ! saving temporary packet data
 IF(MODULO(pack_index, n_pack_save) == 0) THEN
  CALL save_temp_packs(pack_index)
 END IF
END DO
  
END SUBROUTINE update_packages
