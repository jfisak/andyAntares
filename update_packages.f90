 SUBROUTINE update_packages(n_pack, opa_cell, lower_opa, delta_opa)

  USE types

  IMPLICIT NONE    

    INTEGER                           :: I, I_esc, pack_index, nc, next_cell, n_pack
    DOUBLE PRECISION                  :: tau, xi, ran2, tau_rand, cell_dist, event_dist, r,   &
                                         rho_cell, I_beta, opa_cell, lower_opa, delta_opa
    DOUBLE PRECISION, PARAMETER       :: rho=1.D0

!  print*, 'A'
  DO pack_index=1,n_pack
!      print*, 'B'
    IF (MODULO(pack_index,10000) .EQ. 0) print*, 'Working on packet ',pack_index,' ...'
    DO  WHILE (package(pack_index)%active .EQ. 1)
!      print*, 'C'
        IF (package(pack_index)%typ .EQ. type_rpkt) THEN
!      print*, 'D'
!           If the packet is of type rpkt, it represents a photon. So it needs to be propagated.
            CALL do_rpackage(pack_index)
        ELSE IF (package(pack_index)%typ .EQ. type_kpkt) THEN 
!           If the packet is of type kpkt, it represents thermal kinetic energy.
!           Sample all possible cooling processes and randomly select one of them
!            CALL do_kpackage(pack_index)
            print*, 'kpkt found should not happen for now'
        ELSE IF (package(pack_index)%typ .EQ. type_ipkt) THEN 
!           If the packet is of type ipkt, it represents atomic internal energy (excitation/ionization).
!           Calculate all transition probabilities and randomly select one of them (macro-atom formalism)
!           CALL do_ipackage(pack_index)
            print*, 'ipkt found, should not happen for now'

        ELSE
            STOP 'ERROR unknown package typ'
        END IF
    END DO

  END DO

 END SUBROUTINE update_packages
