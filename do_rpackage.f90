SUBROUTINE do_rpackage(pack_index)
 
! Propagation of the photon in 3D grid

USE types
USE constants
USE rates_r

IMPLICIT NONE    

INTEGER                                         :: pack_index, next_cell, event
INTEGER                                         :: get_package_model_index
DOUBLE PRECISION                                :: cell_dist, e_dist

DOUBLE PRECISION, PARAMETER                     :: mininum = 1.E1
INTEGER                                         :: I
DOUBLE PRECISION, DIMENSION(3)                  :: pos, corner, width
DOUBLE PRECISION, DIMENSION(3)                  :: cur_cor, cur_width
DOUBLE PRECISION, DIMENSION(3)                  :: delta_r, cur_pos, cur_corner

INTEGER                                         :: n_thomson
INTEGER                                         :: cur_mgi
! INTEGER                                         :: n_tot_cont
TYPE(rrates)                                    :: actirrates
INTEGER                                         :: pomocna_bunka, cur_pgi
INTEGER                                         :: dummypackage, next_cross
! free free

LOGICAL                                         :: change_of_cell

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
pos = package(pack_index)%pos
cur_cor = dyn_cell(cur_pgi)%corner
cur_width = dyn_cell(cur_pgi)%width

IF(debug == 2) THEN
 CALL find_dyn_cell1(package(pack_index)%pos, pomocna_bunka)
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' cur_pgi = ', cur_pgi, ' neigbors = ', dyn_cell(cur_pgi)%neighbor
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 
 ! write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner/R_sun
 write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner/R_star
!  write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos/R_sun
!  write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_sun
 write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos/R_star
 write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_star

 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
 write(98,*) pos/R_sun, dyn_cell(cur_pgi)%corner/R_sun, dyn_cell(cur_pgi)%width
 ! IF(pos(1) <= cur_cor(1) .or. pos(1) >= cur_cor(1) + cur_width(1) .or. &
 !  pos(2) <= cur_cor(2) .or. pos(2) >= cur_cor(2) + cur_width(2) .or. &
 !  pos(3) <= cur_cor(3) .or. pos(3) >= cur_cor(3) + cur_width(3)) THEN
 !   write(*,*) 'do_rpackage: the packet is out of the propGrid cell'
 !   STOP
 ! END IF
END IF

! if(pack_index == 4) STOP 'do_rpackage: testing'


 CALL boundary3(pack_index, cell_dist, next_cell)
 cur_mgi = get_package_model_index(pack_index)
 next_cross = package(pack_index)%next_cross


 IF((cell_dist > R_inf) .and. (package(pack_index)%virtual .EQV. .FALSE.)) THEN
  write(*,*) 'do_rpackage: cell_dist = ', cell_dist, ' > R_inf'
  write(*,*) 'exiting now'
  STOP
 ELSE IF((cell_dist > R_inf) .and. (package(pack_index)%virtual .EQV. .TRUE.)) THEN
  package(pack_index)%active = 0
 END IF


IF(debug == 2) THEN
 write(*,*) 'do_rpackage 0: cell_dist = ', cell_dist
 ! IF (cell_dist .LT. 0.D0) STOP 'cell_dist < 0'
END IF


IF(cell_dist == 0.e0) STOP 'do_rpackage: cell_dist == 0'

! write(*,*) 'do_rpackage: cur_mgi = ', cur_mgi

IF (cur_mgi > n_modelgrid) THEN
 ! Package is outside the wind model but still inside the propagation grid qube
 ! No physical interaction should occure, set e_dist > cell_dist
 ! write(*,*) 'do_rpackage: vacuum cell for a packet: ', pack_index
 e_dist = cell_dist + R_inf
ELSE
 ! write(*,*) 'do_rpackage: calling event_dist'
 IF(cell_dist > 1.D20) THEN
  write(*,*) 'do_rpackage: pack_index = ', pack_index, ' cell_dist = ', cell_dist
  write(*,*) 'do_rpackage: pos = ', norm2(package(pack_index)%pos)/R_inf
 END IF
 CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
END IF

if(isnan(package(pack_index)%freq_cmf)) STOP 'do_rpackage: freq is NaN'

IF (e_dist .LT. cell_dist) THEN
 change_of_cell = .FALSE.
 CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
 ! if the e_dist == 0 then the packets was re-emmited in the photosphere and no line interaction is allowed
 IF(e_dist > 0.D0) THEN
  CALL update_estimators(pack_index, e_dist)
  CALL do_rpackage_event(pack_index, event, actirrates)
 ! this is an exception when the zero distance is calculated
 ELSE IF (e_dist == 0.0) THEN
  CALL bound_dist(pack_index, cur_pgi, cell_dist)
 END IF
ELSE IF(e_dist > cell_dist .and. cell_dist > 0.e0) THEN 
 ! Move package from the curent position for the cell_dist
 change_of_cell = .TRUE.
 CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
 CALL update_estimators(pack_index, cell_dist)
ELSE IF((cell_dist < 0.e0) .and. (next_cross <= 6) .and. (next_cross >=1)) THEN
 ! next_cell = dyn_cell(cur_pgi)%neighbor(package(pack_index)%next_cross)
 CALL change_cell(pack_index, next_cell)
ELSE IF((cell_dist < 0.e0) .and. (next_cross > 6)) THEN
 cur_corner = dyn_cell(cur_pgi)%corner
 cur_pos = package(pack_index)%pos
 delta_r = (cur_corner - cur_pos)/NORM2(cur_corner - cur_pos)
 if(delta_r(1) == 0.D0) delta_r(1) = 1.D0
 if(delta_r(2) == 0.D0) delta_r(2) = 1.D0
 if(delta_r(3) == 0.D0) delta_r(3) = 1.D0
 IF(next_cross == edyz) THEN
  ! move packet to with the vector (0,1,1)
  package(pack_index)%pos = package(pack_index)%pos + (/0.D0,delta_r(2),delta_r(3)/)
 ELSE IF(next_cross == edxz) THEN
  ! move packet to with the vector (1,0,1)
  package(pack_index)%pos = package(pack_index)%pos + (/delta_r(1),0.D0,delta_r(3)/)
 ELSE IF(next_cross == edxy) THEN
  ! move packet to with the vector (0,1,1)
  package(pack_index)%pos = package(pack_index)%pos + (/delta_r(1),delta_r(2),0.D0/)
 END IF
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
 ! write(*,*) 'do_rpackage II: cell starting = ', dyn_cell(cur_pgi)%corner/R_sun
 ! write(*,*) 'do_rpackage II: packet pos = ', package(pack_index)%pos/R_sun
 ! write(*,*) 'do_rpackage II: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_sun
 write(*,*) 'do_rpackage II: cell starting = ', dyn_cell(cur_pgi)%corner/R_star
 write(*,*) 'do_rpackage II: packet pos = ', package(pack_index)%pos/R_star
 write(*,*) 'do_rpackage II: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/R_star
 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
 DO I = 1,3
  IF(((pos(I) <= corner(I) - mininum) .OR. (pos(I) >= corner(I) + width(I) + mininum)) .and. pack_index /= dummypackage ) THEN
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
