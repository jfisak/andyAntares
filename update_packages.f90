 SUBROUTINE update_packages(n_pack)

  USE types

  IMPLICIT NONE    

  INTEGER             :: n_pack, pack_index
  

  OPEN (UNIT=3, FILE='position.dat')
  ! print*, 'A'
  DO pack_index = 1, n_pack
     !write(*,*) 'update_packages: pack_index = ', pack_index
     IF (MODULO(pack_index,10000) .EQ. 0) print*, 'Working on packet ', pack_index,' ...'
     !print*, 'Working on packet ', pack_index,' ...'
     IF (debug .NE. 0) print*, 'Working on packet ', pack_index,' ...' 
      
     ! Do this loop until something happened with package
     DO  WHILE (package(pack_index)%active .EQ. 1)
        ! print*, 'C'
        IF (package(pack_index)%typ .EQ. type_rpkt) THEN
         IF(package(pack_index)%n_interactions .EQ. 1000000) THEN
          print*, 'package ', pack_index, ' interacted for 2000000 times and will be destroyed...'
          package(pack_index)%active = 0
          destroyed_pack = destroyed_pack + 1
         END IF
           ! print*, 'D'
           ! If the packet is of type rpkt, it represents a photon. So it needs to be propagated.
           CALL do_rpackage(pack_index)
        ELSE IF (package(pack_index)%typ .EQ. type_kpkt) THEN 
           ! If the packet is of type kpkt, it represents thermal kinetic energy.
           ! Sample all possible cooling processes and randomly select one of them
            CALL do_kpackage(pack_index)
           !print*, 'kpkt found should not happen for now'
        ELSE IF (package(pack_index)%typ .EQ. type_ipkt) THEN 
           ! If the packet is of type ipkt, it represents atomic internal energy (excitation/ionization).
           ! Calculate all transition probabilities and randomly select one of them (macro-atom formalism)
           CALL do_ipackage(pack_index)
           ! print*, 'ipkt found, should not happen for now'
        ELSE
           STOP 'ERROR unknown package typ'
        END IF
     END DO

  END DO
  CLOSE(UNIT=3)
  
END SUBROUTINE update_packages
