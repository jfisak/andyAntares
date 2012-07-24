SUBROUTINE do_rpackage(pack_index)
 
  ! Propagation of the photon in 3D grid

  USE types

  IMPLICIT NONE    

  INTEGER                           :: I, I_esc, pack_index, nc, next_cell, n_pack, event
  DOUBLE PRECISION                  :: tau, xi, ran2, tau_rand, cell_dist, e_dist, r,   &
                                         rho_cell, I_beta, opa_cell, lower_opa, delta_opa
  DOUBLE PRECISION, PARAMETER       :: rho=1.D0

  ! DOUBLE PRECISION, PARAMETER       :: opa_cell=1.D0/20.D0, rho=1.D0
  ! DOUBLE PRECISION, PARAMETER       :: opa_cell=5.D-3, rho=5.D-2

  OPEN (UNIT=3, FILE='position.dat')  
  ! WRITE(3,*) pack_index, package(pack_index)%cell_numb, package(pack_index)%pos


  CALL boundary(pack_index, cell_dist, next_cell)
  IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
  CALL event_dist(pack_index, cell_dist, e_dist, event)
  IF (debug .EQ. 1) THEN 
     print*, cell_dist, next_cell, e_dist , cell(package(pack_index)%cell_numb)%indexc
  END IF



  IF (e_dist .LT. cell_dist) THEN   
     ! Move photon package from the curent position for some distance
     !print*, pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
     CALL move_package(pack_index, e_dist)
     !print*, pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
     CALL do_rpackage_event(pack_index, event)
     IF (debug .EQ. 1) THEN 
        print*, 'do event', opa_cell * rho_cell * cell_dist
     END IF
  ELSE     
     ! Move package from the curent position for the cell_dist
     CALL move_package(pack_index, cell_dist)
     ! If package escaped the calculation volume (next_cell=-99) then
     ! it become no-active and package type is update to the
     ! type_escaped, else the cell number is updated
     CALL change_cell(pack_index, next_cell)
     IF (debug .EQ. 1) THEN 
        print*, 'propagate ', opa_cell * rho_cell * cell_dist
     END IF
  END IF
  

END SUBROUTINE do_rpackage
