SUBROUTINE do_rpackage(pack_index)
 
  ! Propagation of the photon in 3D grid

  USE types
  USE rates_r

  IMPLICIT NONE    

  INTEGER                           :: pack_index, next_cell, event
  INTEGER                           :: get_package_model_index
  DOUBLE PRECISION                  :: cell_dist, e_dist, &
                                       rho_cell

INTEGER                                         :: n_thomson
INTEGER                                         :: cur_mgi
!INTEGER                                         :: n_tot_cont
TYPE(rrates)                                    :: actirrates
! free free

LOGICAL                                 :: change_of_cell

! number of thomson scattering 
! total number of continuum opacity sources
IF(n_ff /= 0) THEN
! nothing to be done in this fork
ELSE
  n_ff = 1
  n_thomson = 1
  n_tot_cont = n_thomson + n_photcrossect + n_ff
END IF
actirrates = rrates()

  CALL boundary3(pack_index, cell_dist, next_cell)
  IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
  cur_mgi = get_package_model_index(pack_index)
  IF (cur_mgi .EQ. n_modelgrid + 1) THEN
   ! Package is outside the wind model but still inside the propagation grid qube
   ! No physical interaction should occure, set e_dist > cell_dist
   e_dist = cell_dist + R_inf
  ELSE IF(cur_mgi .EQ. n_modelgrid + 2) THEN
   e_dist = 1.D50
  ELSE
   ! write(*,*) 'do_rpackage: calling event_dist'
   IF(cell_dist > 1.D20) write(*,*) 'do_rpackage: pack_index = ', pack_index, ' cell_dist = ', cell_dist
   CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
  END IF

  IF (e_dist .LT. cell_dist) THEN
   change_of_cell = .FALSE.
   CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
   ! if the e_dist == 0 then the packets was re-emmited in the photosphere and no line interaction is allowed
   IF(e_dist > 0.D0) THEN
    CALL update_estimators(pack_index, e_dist)
    CALL do_rpackage_event(pack_index, event, actirrates)
   ELSE IF (e_dist == 0.0) THEN
    CALL find_dist(pack_index, next_cell, cell_dist)
   END IF
  ELSE     
   ! Move package from the curent position for the cell_dist
   change_of_cell = .TRUE.
   CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
   CALL update_estimators(pack_index, cell_dist)
  END IF
  
END SUBROUTINE do_rpackage
