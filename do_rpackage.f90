SUBROUTINE do_rpackage(pack_index)
 
! Propagation of the photon in 3D grid

USE types
USE constants
USE rates_r
USE counters

IMPLICIT NONE    

INTEGER                                         :: pack_index, next_cell, event
INTEGER                                         :: get_package_model_index
DOUBLE PRECISION                                :: cell_dist, e_dist

DOUBLE PRECISION, PARAMETER                     :: mininum = 1.E-1
INTEGER                                         :: ind_I
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos, corner, width
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: cur_cor, cur_width
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: delta_r, cur_pos, cur_corner

INTEGER                                         :: n_thomson
INTEGER                                         :: cur_mgi
! INTEGER                                         :: n_tot_cont
TYPE(rrates)                                    :: actirrates
INTEGER                                         :: pomocna_bunka, cur_pgi
INTEGER                                         :: next_cross
! DOUBLE PRECISION                                :: max_dist
LOGICAL, PARAMETER                              :: procout = .TRUE.
! free free

LOGICAL                                         :: change_of_cell

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


! initial variables
cur_pgi = package(pack_index)%cell_numb
pos = package(pack_index)%pos
cur_cor = dyn_cell(cur_pgi)%corner
cur_width = dyn_cell(cur_pgi)%width
if(pack_index == 21) write(98,*) pos/const_Rsun, dyn_cell(cur_pgi)%corner/const_Rsun, dyn_cell(cur_pgi)%width/const_Rsun

! the first debug part
IF(debug == 2) THEN
 CALL find_dyn_cell1(package(pack_index)%pos, pomocna_bunka)
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' cur_pgi = ', cur_pgi, ' neigbors = ', dyn_cell(cur_pgi)%neighbor
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 
 ! write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner/const_Rsun
 write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner
!  write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos/const_Rsun
!  write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)/const_Rsun
 write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos
 write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)

 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
 cur_pgi = package(pack_index)%cell_numb
 pos = package(pack_index)%pos
 corner = dyn_cell(cur_pgi)%corner
 width = dyn_cell(cur_pgi)%width
 CALL find_dyn_cell1(pos, pomocna_bunka)
 write(*,*) 'do_rpackage I: cur_pgi = ', cur_pgi
 write(*,*) 'do_rpackage I: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 write(*,*) 'do_rpackage I: cur_pgi = ', cur_pgi
 write(*,*) 'do_rpackage I: cell starting = ', dyn_cell(cur_pgi)%corner
 write(*,*) 'do_rpackage I: packet pos = ', package(pack_index)%pos
 write(*,*) 'do_rpackage I: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)
 write(*,*) 'do_rpackage I: direction = ', package(pack_index)%dir
 write(*,*) 'do_rpackage I: active = ', package(pack_index)%active
END IF


! looking for the next boundary which the packet crosses as the first and an index of this boundary
CALL boundary3(pack_index, cell_dist, next_cell)
cur_mgi = get_package_model_index(pack_index)
next_cross = package(pack_index)%next_cross
write(*,*) 'do_rpackage: cell_dist = ', cell_dist
IF(procout) write(*,*) 'event_dist: after calling boundary3, cell_dist = ', cell_dist, ' next_cell = ', next_cell


! testing if cell_dist is still inside the propGrid
IF((cell_dist > R_inf) .and. (package(pack_index)%virtual .EQV. .FALSE.)) THEN
 write(*,*) 'do_rpackage: cell_dist = ', cell_dist, ' > R_inf'
 write(*,*) 'exiting now'
 STOP
ELSE IF((cell_dist > R_inf) .and. (package(pack_index)%virtual .EQV. .TRUE.)) THEN
 package(pack_index)%active = 0
END IF

IF(cell_dist == 0.e0) STOP 'do_rpackage: cell_dist == 0'

! test on cur_mgi variable which should be <= n_modelgrid, othervise no
! interaction is allowed
IF(cur_mgi <= n_modelgrid) THEN
 IF(procout) write(*,*) 'do_rpackage: calling event_dist'
 CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
 IF(procout) write(*,*) 'do_rpackage: e_dist/bdist = ', e_dist/cell_dist, ' event = ', event
ELSE IF(cur_mgi > n_modelgrid) THEN
 ! Package is outside the wind model but still inside the propagation grid cube
 ! No physical interaction should occure, set e_dist > cell_dist
 ! write(*,*) 'do_rpackage: vacuum cell for a packet: ', pack_index
 IF(procout) write(*,*) 'event_dist: packet is outside the modGrid, cur_mgi = ', cur_mgi, n_modelgrid
 e_dist = cell_dist + R_inf
END IF

if(isnan(package(pack_index)%freq_cmf)) STOP 'do_rpackage: freq is NaN'


IF (e_dist < cell_dist) THEN
 change_of_cell = .FALSE.
 IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
 CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
 ! if the e_dist == 0 then the packets was re-emmited in the photosphere and no line interaction is allowed
 IF(e_dist > 0.D0) THEN
  IF(procout) write(*,*) 'do_rpackage: move_package, e_dist = ', e_dist, ' > 0.D0'
  CALL update_estimators(pack_index, e_dist)
  CALL do_rpackage_event(pack_index, event, actirrates)
 ! this is an exception when the zero distance is calculated
 ELSE IF (e_dist == 0.0) THEN
  IF(procout) write(*,*) 'do_rpackage: move_package, e_dist = 0.D0'
  CALL bound_dist(pack_index, cur_pgi, cell_dist)
 END IF
ELSE IF(e_dist > cell_dist .and. cell_dist > 0.e0) THEN 
 ! Move package from the curent position for the cell_dist
 change_of_cell = .TRUE.
 IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
 CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
 CALL update_estimators(pack_index, cell_dist)
ELSE IF((cell_dist < 0.e0) .and. (next_cross <= 6) .and. (next_cross >=1)) THEN
 ! next_cell = dyn_cell(cur_pgi)%neighbor(package(pack_index)%next_cross)
 CALL change_cell(pack_index, next_cell)
ELSE IF((cell_dist < 0.e0) .and. (next_cross > 6)) THEN
 cur_corner = dyn_cell(cur_pgi)%corner
 cur_pos = package(pack_index)%pos
 delta_r = (cur_corner - cur_pos)/NORM2(cur_corner - cur_pos)
 if(delta_r(1) == 0.D0) delta_r(1) = 1.D1
 if(delta_r(2) == 0.D0) delta_r(2) = 1.D1
 if(delta_r(3) == 0.D0) delta_r(3) = 1.D1
 IF(next_cross == edyz) THEN
  ! move packet to with the vector (0,1,1)
  package(pack_index)%pos = package(pack_index)%pos + (/0.D0,delta_r(2),delta_r(3)/)
  IF(procout) write(*,*) 'do_rpackage: move packet with the vector (0,1,1)'
 ELSE IF(next_cross == edxz) THEN
  ! move packet to with the vector (1,0,1)
  package(pack_index)%pos = package(pack_index)%pos + (/delta_r(1),0.D0,delta_r(3)/)
  IF(procout) write(*,*) 'do_rpackage: move packet with the vector (1,0,1)'
 ELSE IF(next_cross == edxy) THEN
  ! move packet to with the vector (0,1,1)
  package(pack_index)%pos = package(pack_index)%pos + (/delta_r(1),delta_r(2),0.D0/)
  IF(procout) write(*,*) 'do_rpackage: move packet with the vector (0,1,1)'
 END IF
END IF

cur_pgi = package(pack_index)%cell_numb
pos = package(pack_index)%pos
corner = dyn_cell(cur_pgi)%corner
width = dyn_cell(cur_pgi)%width
write(*,*) 'do_rpackage III: cell starting = ', dyn_cell(cur_pgi)%corner
write(*,*) 'do_rpackage III: packet pos = ', pos
write(*,*) 'do_rpackage III: cell ending = ', (corner + width)
DO ind_I = 1, const_dimofspace
 IF(((pos(ind_I) < corner(ind_I)) .OR. (pos(ind_I) > corner(ind_I) + width(ind_I))) &
  .and. pack_index <= SIZE(package) ) THEN
  if(pack_index == 21) write(98,*) pos/const_Rsun, dyn_cell(cur_pgi)%corner/const_Rsun, dyn_cell(cur_pgi)%width/const_Rsun
  CALL correction_propagation(ind_I, pack_index, cur_pgi)
  if(pack_index == 21) write(98,*) pos/const_Rsun, dyn_cell(cur_pgi)%corner/const_Rsun, dyn_cell(cur_pgi)%width/const_Rsun
  cur_pgi = package(pack_index)%cell_numb
  pos = package(pack_index)%pos
  corner = dyn_cell(cur_pgi)%corner
  width = dyn_cell(cur_pgi)%width
 END IF
END DO
write(*,*) 'do_rpackage III: cell starting = ', dyn_cell(cur_pgi)%corner
write(*,*) 'do_rpackage III: packet pos = ', package(pack_index)%pos
write(*,*) 'do_rpackage III: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)
! check if we are in a correct cell
IF(debug == 2) THEN
 cur_pgi = package(pack_index)%cell_numb
 pos = package(pack_index)%pos
 corner = dyn_cell(cur_pgi)%corner
 width = dyn_cell(cur_pgi)%width
 write(*,*) 'do_rpackage: pos = ', pos, ' corner = ', corner, ' width = ', width
 CALL find_dyn_cell1(pos, pomocna_bunka)
 write(*,*) 'do_rpackage II: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 write(*,*) 'do_rpackage II: cur_pgi = ', cur_pgi
 write(*,*) 'do_rpackage II: cell starting = ', dyn_cell(cur_pgi)%corner
 write(*,*) 'do_rpackage II: packet pos = ', package(pack_index)%pos
 write(*,*) 'do_rpackage II: cell ending = ', (dyn_cell(cur_pgi)%corner + dyn_cell(cur_pgi)%width)
 write(*,*) 'do_rpackage II: direction = ', package(pack_index)%dir
 write(*,*) 'do_rpackage II: active = ', package(pack_index)%active
 DO ind_I = 1, const_dimofspace
  IF(((pos(ind_I) < corner(ind_I) - mininum) .OR. (pos(ind_I) > corner(ind_I) + width(ind_I) + mininum)) &
   .and. pack_index <= SIZE(package)) THEN
   write(*,*) 'do_rpackage i: r - min = ', (corner(ind_I) - mininum), ' r + w + min = ', (corner(ind_I) + width(ind_I) + mininum)
   CALL find_dyn_cell1(pos, pomocna_bunka)
   ! write(*,*) 'do_rpackage: pos - cell = ', 
   write(*,*) 'do_rpackage II: cell starting = ', corner
   write(*,*) 'do_rpackage II: packet pos = ', pos
   write(*,*) 'do_rpackage II: cell ending = ', corner + width
   write(*,*) 'do_rpackage: skutecna bunka = ', pomocna_bunka
   write(*,*) 'do_rpackage ii: r - min = ', (corner(ind_I) - mininum), ' r + w + min = ', (corner(ind_I) + width(ind_I) + mininum)
   write(4,*) pos, dyn_cell(cur_pgi)%corner, dyn_cell(cur_pgi)%width, get_package_model_index(pack_index)
   write(4,*) pos, corner, width, get_package_model_index(pack_index)
   write(*,*) 'I = ', ind_I
   STOP 'do_rpackage: packet is not located inside the propagation cell'
  END IF
 END DO
END IF

  
END SUBROUTINE do_rpackage
