SUBROUTINE do_rpackage(pack_index)
 
  ! Propagation of the photon in 3D grid

  USE types
  USE rates_r

  IMPLICIT NONE    

  INTEGER                           :: pack_index, next_cell, event
  INTEGER                           :: get_package_model_index
  DOUBLE PRECISION                  :: cell_dist, e_dist, &
                                       rho_cell, opa_cell
  ! DOUBLE PRECISION, PARAMETER       :: rho = 1.D0

  ! DOUBLE PRECISION, PARAMETER       :: opa_cell=1.D0/20.D0, rho=1.D0
  ! DOUBLE PRECISION, PARAMETER       :: opa_cell=5.D-3, rho=5.D-2
INTEGER                                         :: n_thomson
!INTEGER                                         :: n_tot_cont
TYPE(rrates)                                    :: actirrates
! free free

! number of thomson scattering 
! total number of continuum opacity sources
!write(*,*) 'do_rpackage: n_photcrossect = ', n_photcrossect
IF(n_ff /= 0) THEN
! nothing to be done in this fork
ELSE
  n_ff = 1
  n_thomson = 1
  n_tot_cont = n_thomson + n_photcrossect + n_ff
END IF
actirrates = rrates()

  CALL boundary3(pack_index, cell_dist, next_cell)
  ! WRITE(3,*) pack_index, package(pack_index)%pos/R_sun, dyn_cell(package(pack_index)%cell_numb)%corner/R_sun, &
  !               dyn_cell(package(pack_index)%cell_numb)%width/R_sun, get_package_model_index(pack_index)
  ! write(78,*) pack_index, package(pack_index)%freq_rf, package(pack_index)%freq_cmf, package(pack_index)%freq_cmf / linelist(1)%freq
  IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
  IF (get_package_model_index(pack_index) .EQ. n_modelgrid + 1) THEN
      ! Package is outside the wind model but still inside the propagation grid qube
      ! No physical interaction should occure, set e_dist > cell_dist
      e_dist = cell_dist + 1.D10
  ELSE IF(get_package_model_index(pack_index) .EQ. n_modelgrid + 2) THEN
      e_dist = 1.D50
      !write(99,*) 'package: ', pack_index, ' is in empty space...'
  ELSE
      ! write(*,*) 'do_rpackage: calling event_dist'
      CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
  END IF

  ! write(99,*) 'e_dist = ', e_dist, ' cell_dist = ', cell_dist
  IF (e_dist .LT. cell_dist) THEN
     ! write(99,*) 'photon interacts'   
     ! Move photon package from the current position for some distance
     ! write(99,*) pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
     !write(99,*) 'before moving package #', pack_index
     ! write(*,*) 'before moving package #', pack_index
     CALL move_package(pack_index, e_dist)
     CALL update_estimators(pack_index, e_dist)
     ! write(99,*) pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
     CALL do_rpackage_event(pack_index, event, actirrates)
     IF (debug .EQ. 1) THEN 
        write(99,*) 'do event', opa_cell * rho_cell * cell_dist
     END IF
  ELSE     
     ! Move package from the curent position for the cell_dist
     ! write(*,*) 'before moving package ##', pack_index
     !write(99,*) 'before moving package #', pack_index
     CALL move_package(pack_index, cell_dist)     
     CALL update_estimators(pack_index, cell_dist)
     ! If package escaped the calculation volume (next_cell=-99) then
     ! it become no-active and package type is update to the
     ! type_escaped, else the cell number is updated
     CALL change_cell(pack_index, next_cell)
     IF (debug .EQ. 1) THEN 
        write(99,*) 'propagate ', opa_cell * rho_cell * cell_dist
     END IF
  END IF
  
END SUBROUTINE do_rpackage
