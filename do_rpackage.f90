SUBROUTINE do_rpackage(pack_index)
 
! Propagation of the photon in 3D grid

USE types
USE rates_r

IMPLICIT NONE    

INTEGER                           :: pack_index, next_cell, event
INTEGER                           :: get_package_model_index
DOUBLE PRECISION                  :: cell_dist, e_dist, &
                                     rho_cell

DOUBLE PRECISION, PARAMETER           :: mininum = 1.E1
INTEGER                                 :: I
DOUBLE PRECISION, DIMENSION(3)          :: pos, corner, width

INTEGER                                         :: n_thomson
INTEGER                                         :: cur_mgi
!INTEGER                                         :: n_tot_cont
TYPE(rrates)                                    :: actirrates
INTEGER                                         :: pomocna_bunka, cur_pgi
INTEGER                                 :: dummypackage
! free free

LOGICAL                                 :: change_of_cell

dummypackage = SIZE(package)

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

cur_pgi = package(pack_index)%cell_numb

! write(23,*) package(pack_index)%pos

IF(debug == 2) THEN
 CALL find_dyn_cell1(package(pack_index)%pos, pomocna_bunka)
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' cur_pgi = ', cur_pgi, ' neigbors = ', dyn_cell(cur_pgi)%neighbor
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 
 write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner/R_sun
 write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos/R_sun
 write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_sun

 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
END IF



 CALL boundary3(pack_index, cell_dist, next_cell)
 cur_mgi = get_package_model_index(pack_index)



IF(debug == 2) THEN
 write(*,*) 'do_rpackage 0: cell_dist = ', cell_dist
 ! IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
END IF


IF(cell_dist == 0.e0) STOP 'do_rpackage: cell_dist == 0'


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
  CALL find_dist(pack_index, cur_pgi, cell_dist)
 END IF
ELSE IF(e_dist > cell_dist .and. cell_dist > 0.e0) THEN 
 ! Move package from the curent position for the cell_dist
 change_of_cell = .TRUE.
 CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
 CALL update_estimators(pack_index, cell_dist)
ELSE IF(cell_dist < 0.e0) THEN
 ! next_cell = dyn_cell(cur_pgi)%neighbor(package(pack_index)%next_cross)
 CALL change_cell(pack_index, next_cell)
END IF

! check if we are in a correct cell
IF(debug == 2) THEN
 CALL find_dyn_cell1(package(pack_index)%pos, pomocna_bunka)
 cur_pgi = package(pack_index)%cell_numb
 pos = package(pack_index)%pos
 corner = dyn_cell(cur_pgi)%corner
 width = dyn_cell(cur_pgi)%width
 write(*,*) 'do_rpackage II: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 write(*,*) 'do_rpackage II: cur_pgi = ', cur_pgi
 write(*,*) 'do_rpackage II: cell starting = ', dyn_cell(cur_pgi)%corner/R_sun
 write(*,*) 'do_rpackage II: packet pos = ', package(pack_index)%pos/R_sun
 write(*,*) 'do_rpackage II: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_sun
 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
 DO I = 1,3
  IF((pos(I) < corner(I) - mininum .OR. pos(I) > corner(I) + width(I) + mininum) .and. pack_index /= dummypackage ) THEN
   CALL find_dyn_cell1(pos, pomocna_bunka)
   write(*,*) 'do_rpackage: skutecna bunka = ', pomocna_bunka
   write(4,*) pos, dyn_cell(cur_pgi)%corner, dyn_cell(cur_pgi)%width, get_package_model_index(pack_index)
   write(4,*) pos, corner, width, get_package_model_index(pack_index)
   write(*,*) 'I = ', I
   STOP 'do_rpackage: packet is not located inside the propagation cell'
  END IF
 END DO
END IF

  
END SUBROUTINE do_rpackage
