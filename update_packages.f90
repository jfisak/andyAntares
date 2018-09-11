 SUBROUTINE update_packages(n_pack)

  USE types
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
 IF (MODULO(pack_index,100000) .EQ. 0) write(99,*) 'Working on packet ', pack_index,' ...'
 ! IF (MODULO(pack_index,100000) .EQ. 0) write(*,*) 'Working on packet ', pack_index,' ...'
 ! write(*,*) 'Working on packet ', pack_index,' ...'
 ! write(*,*) 'Working on packet ', pack_index,' ...'
 ! Do this loop until something happened with package
 DO  WHILE (package(pack_index)%active .EQ. 1)
  ! write(*,*) 'C'
  IF (package(pack_index)%typ .EQ. type_rpkt) THEN
   IF(package(pack_index)%n_interactions .EQ. 1000000) THEN
    ! write(*,*) 'package ', pack_index, ' interacted for 2000000 times and will be destroyed...'
    package(pack_index)%active = 0
    !$OMP ATOMIC
    count_des_inte = count_des_inte + 1
   END IF
     ! write(*,*) 'D'
     ! If the packet is of type rpkt, it represents a photon. So it needs to be propagated.
     ! write(*,*) 'update_packages: calling do_rpackage'
     CALL do_rpackage(pack_index)
  ELSE IF (package(pack_index)%typ .EQ. type_kpkt) THEN 
     ! If the packet is of type kpkt, it represents thermal kinetic energy.
     ! Sample all possible cooling processes and randomly select one of them
     ! write(*,*) 'update_packages: calling do_kpackage'
      CALL do_kpackage(pack_index)
     !write(*,*) 'kpkt found should not happen for now'
  ELSE IF (package(pack_index)%typ .EQ. type_ipkt) THEN 
     ! If the packet is of type ipkt, it represents atomic internal energy (excitation/ionization).
     ! Calculate all transition probabilities and randomly select one of them (macro-atom formalism)
     ! write(*,*) 'update_packages: calling do_ipackage'
     CALL do_ipackage(pack_index)
     ! write(*,*) 'ipkt found, should not happen for now'
  ELSE
     STOP 'ERROR unknown package typ'
  END IF
 END DO
 ! saving temporary packet data
 IF(MODULO(pack_index, n_pack_save) == 0) THEN
  CALL save_temp_packs(pack_index)
 END IF
END DO
  
END SUBROUTINE update_packages
