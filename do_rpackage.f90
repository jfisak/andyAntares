SUBROUTINE do_rpackage(pack_index)
 
  ! Propagation of the photon in 3D grid

  USE types
  USE rates_r

  IMPLICIT NONE    

  INTEGER                           :: I, I_esc, pack_index, nc, next_cell, n_pack, event
  INTEGER                           :: get_package_model_index
  DOUBLE PRECISION                  :: tau, xi, ran2, tau_rand, cell_dist, e_dist, r,   &
                                       rho_cell, I_beta, opa_cell, lower_opa, delta_opa
  INTEGER                           :: check_cell
  ! DOUBLE PRECISION, PARAMETER       :: rho = 1.D0

  ! DOUBLE PRECISION, PARAMETER       :: opa_cell=1.D0/20.D0, rho=1.D0
  ! DOUBLE PRECISION, PARAMETER       :: opa_cell=5.D-3, rho=5.D-2
INTEGER                                         :: n_thomson
!INTEGER                                         :: n_tot_cont
INTEGER                                         :: my_rank
TYPE(rrates)                                    :: actirrates

actirrates = rrates()
! number of thomson scattering (code must be clear)
n_thomson = 1
! total number of continuum opacity sources
n_tot_cont = n_thomson + n_photcrossect

  CALL boundary3(pack_index, cell_dist, next_cell)
!  WRITE(3,*) package(pack_index)%pos, dyn_cell(package(pack_index)%cell_numb)%corner, &
!                dyn_cell(package(pack_index)%cell_numb)%width, get_package_model_index(pack_index)
  IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
  IF (get_package_model_index(pack_index) .EQ. n_modelgrid + 1) THEN
      ! Package is outside the wind model but still inside the propagation grid qube
      ! No physical interaction should occure, set e_dist > cell_dist
      e_dist = cell_dist + 1.D10
  ELSE IF(get_package_model_index(pack_index) .EQ. n_modelgrid + 2) THEN
      e_dist = 1.D50
      !print*, 'package: ', pack_index, ' is in empty space...'
  ELSE
      CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
  END IF

!  IF (debug .EQ. 1) THEN 
!      print*, cell_dist, next_cell, e_dist , cell(package(pack_index)%cell_numb)%model_index
!  END IF

  !print*, 'e_dist = ', e_dist, ' cell_dist = ', cell_dist
  IF (e_dist .LT. cell_dist) THEN
     ! print*, 'photon interacts'   
     ! Move photon package from the current position for some distance
     ! print*, pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
     !print*, 'before moving package #', pack_index
     CALL move_package(pack_index, e_dist)
     CALL update_estimators(pack_index, e_dist)
     ! print*, pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
     CALL do_rpackage_event(pack_index, event, actirrates)
     IF (debug .EQ. 1) THEN 
        print*, 'do event', opa_cell * rho_cell * cell_dist
     END IF
  ELSE     
     ! Move package from the curent position for the cell_dist
     !print*, 'before moving package #', pack_index
     CALL move_package(pack_index, cell_dist)     
     CALL update_estimators(pack_index, cell_dist)
     ! If package escaped the calculation volume (next_cell=-99) then
     ! it become no-active and package type is update to the
     ! type_escaped, else the cell number is updated
     CALL change_cell(pack_index, next_cell)
     IF (debug .EQ. 1) THEN 
        print*, 'propagate ', opa_cell * rho_cell * cell_dist
     END IF
  END IF
  

END SUBROUTINE do_rpackage
