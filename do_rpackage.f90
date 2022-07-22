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
  ! write(*,*) 'do_rpackage: cur_index = ', package(pack_index)%cell_numb, ' next_cell = ', next_cell
   ! IF(pack_index == 78) THEN
   !  WRITE(3,*) package(pack_index)%pos, dyn_cell(package(pack_index)%cell_numb)%corner, &
   !                dyn_cell(package(pack_index)%cell_numb)%width, get_package_model_index(pack_index)
   ! END IF
  IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
  cur_mgi = get_package_model_index(pack_index)
  IF (cur_mgi .EQ. n_modelgrid + 1) THEN
      ! Package is outside the wind model but still inside the propagation grid qube
      ! No physical interaction should occure, set e_dist > cell_dist
      e_dist = cell_dist + R_inf
  ELSE IF(cur_mgi .EQ. n_modelgrid + 2) THEN
      e_dist = 1.D50
      !write(99,*) 'package: ', pack_index, ' is in empty space...'
  ELSE
      ! write(*,*) 'do_rpackage: calling event_dist'
      IF(cell_dist > 1.D20) write(*,*) 'do_rpackage: pack_index = ', pack_index, ' cell_dist = ', cell_dist
      CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
  END IF
  ! write(*,*) 'do_rpackage: e_dist = ', e_dist/R_sun, ' cell_dist = ', cell_dist/R_sun, ' cur_mgi = ', cur_mgi, n_modelgrid

  ! write(99,*) 'e_dist = ', e_dist, ' cell_dist = ', cell_dist
  IF (e_dist .LT. cell_dist) THEN
   ! write(99,*) 'photon interacts'   
   ! Move photon package from the current position for some distance
   ! write(99,*) pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
   !write(99,*) 'before moving package #', pack_index
   ! write(*,*) 'before moving package #', pack_index
   ! write(*,*) 'do_rpackage: #1 cur_index = ', package(pack_index)%cell_numb, ' next_cell = ', next_cell
   ! write(*,*) 'before moving package ##', pack_index
   change_of_cell = .FALSE.
   CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
   ! write(*,*) 'after moving package ##', pack_index
   IF(e_dist > 0.D0) THEN
    CALL update_estimators(pack_index, e_dist)
    ! write(99,*) pack_index, freq_line, package(pack_index)%freq_cmf, package(pack_index)%freq_rf
    CALL do_rpackage_event(pack_index, event, actirrates)
   END IF
  ELSE     
     ! Move package from the curent position for the cell_dist
     !write(99,*) 'before moving package #', pack_index
     ! write(*,*) 'do_rpackage: #2 cur_index = ', package(pack_index)%cell_numb, ' next_cell = ', next_cell
     ! write(*,*) 'before moving package ##', pack_index
     change_of_cell = .TRUE.
     CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
     ! write(*,*) 'after moving package ##', pack_index
     CALL update_estimators(pack_index, cell_dist)
     ! If package escaped the calculation volume (next_cell=-99) then
     ! it become no-active and package type is update to the
     ! type_escaped, else the cell number is updated
     ! write(*,*) 'do_rpackage: e_dist = ', cell_dist
     ! IF(cell_dist > 0.D0) THEN
     !  CALL change_cell(pack_index, next_cell)
     ! END IF
  END IF
  
END SUBROUTINE do_rpackage
