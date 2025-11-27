! Propagation of the r-package in the propGrid
!
! INPUT: pack_index(INT): index of a packet
! OUTPUT: NONE
!
! 1x RETURN POINT
! 
SUBROUTINE do_rpackage(pack_index)
 

USE types
USE constants
USE rates_r

IMPLICIT NONE    

INTEGER                                         :: pack_index, next_cell, event
INTEGER                                         :: get_package_model_index
DOUBLE PRECISION                                :: cell_dist, e_dist

DOUBLE PRECISION, PARAMETER                     :: mininum = 5.E0
INTEGER                                         :: ind_I
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos, corner, upcorner
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: cur_cor, cur_width
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: delta_r, cur_pos, cur_corner, cur_upcorner

INTEGER                                         :: n_thomson
INTEGER                                         :: cur_mgi
! INTEGER                                         :: n_tot_cont
TYPE(rrates)                                    :: actirrates
INTEGER                                         :: pomocna_bunka, cur_pgi
INTEGER                                         :: next_cross
! DOUBLE PRECISION                                :: max_dist
LOGICAL                                         :: procout = .FALSE.
DOUBLE PRECISION, PARAMETER                     :: epsilon0 = 1.D0
! free free
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: pack_pos

LOGICAL                                         :: change_of_cell

INTEGER                                         :: n_pos, n_neg, n_zer, n_par

IF(debug == 2) procout = .TRUE.
!______________________________________________________________________________
! setting up the counters
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


!______________________________________________________________________________
! initial position and current propGrid 
cur_pgi = package(pack_index)%cell_numb
cur_mgi = get_package_model_index(pack_index)
pos = package(pack_index)%pos
cur_cor = dyn_cell(cur_pgi)%corner
cur_width = dyn_cell(cur_pgi)%upcorner - dyn_cell(cur_pgi)%corner

! write(*,*) 'do_rpackage: ||pos||/R_star = ', norm2(pos)/R_star, R_inf, R_star

!______________________________________________________________________________
! seeking for the next boundary (propGrid) and 
! looking for the next boundary which the packet crosses as the first and an index of this boundary
CALL boundary3(pack_index, cell_dist, next_cell, n_pos, n_neg, n_zer, n_par)
next_cross = package(pack_index)%next_cross
IF(procout) write(*,*) 'do_rpackage: after calling boundary3, cell_dist = ', cell_dist, ' next_cell = ', next_cell


!______________________________________________________________________________
! testing if cell_dist is still inside the propGrid
IF((cell_dist > R_inf) .and. (package(pack_index)%virtual .EQV. .FALSE.)) THEN
 write(*,*) 'do_rpackage: cell_dist = ', cell_dist, ' > R_inf'
 write(*,*) 'exiting now'
 package(pack_index)%active = 0
 ! RETURN POINT
 RETURN
 STOP
ELSE IF((cell_dist > R_inf) .and. (package(pack_index)%virtual .EQV. .TRUE.)) THEN
 package(pack_index)%active = 0
END IF

! IF(cell_dist == 0.e0) STOP 'do_rpackage: cell_dist == 0'

!______________________________________________________________________________
! test on cur_mgi variable which should be <= n_modelgrid, othervise no
! interaction is allowed
IF(cur_mgi <= n_modelgrid) THEN
 IF(procout) write(*,*) 'do_rpackage: calling event_dist'
 IF(procout) write(*,*) 'do_rpackage: is diff = ', model_grid(cur_mgi)%is_difapp
 IF(procout) write(*,*) 'do_rpackage: pack type = ', package(pack_index)%typ
 CALL event_dist(pack_index, cell_dist, e_dist, event, actirrates)
 IF(procout) write(*,*) 'do_rpackage: e_dist = ', e_dist
 IF(procout) write(*,*) 'do_rpackage: e_dist/bdist = ', e_dist/cell_dist, ' event = ', event
 ! IF(e_dist < 1.D0) THEN
 !  write(50,*) pos, package(pack_index)%dir, corner
 ! END IF
ELSE IF(cur_mgi > n_modelgrid) THEN
 ! Package is outside the wind model but still inside the propagation grid cube
 ! No physical interaction should occure, set e_dist > cell_dist
 ! write(*,*) 'do_rpackage: vacuum cell for a packet: ', pack_index
 IF(procout) write(*,*) 'event_dist: packet is outside the modGrid, cur_mgi = ', cur_mgi, n_modelgrid
 e_dist = cell_dist + R_inf
END IF

if(isnan(package(pack_index)%freq_cmf)) STOP 'do_rpackage: freq is NaN'


!_____________________________________________________________________________________________________
! now decide what happens with the packet
!_____________________________________________________________________________________________________
IF(procout) write(*,*) 'do_rpackage: e_dist = ', e_dist, ' cell_dist = ', cell_dist, &
 ' n_pos = ', n_pos, ' n_neg = ', n_neg, ' n_zer = ', n_zer, ' n_par = ', n_par
IF (e_dist < cell_dist) THEN
 IF(procout) write(*,*) 'do_rpackage: e_dist < cell_dist'
 IF(n_pos == 3 .and. n_neg == 2 .and. n_zer == 1) THEN
  change_of_cell = .FALSE.
  IF(procout) write(*,*) 'do_rpackage: interaction, n_pos = 3, n_neg = 2, n_zer = 1'
  IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
  CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
 ELSE IF(n_pos == 2 .and. n_neg == 3 .and. n_zer == 1) THEN
  change_of_cell = .FALSE.
  IF(procout) write(*,*) 'do_rpackage: interaction, n_pos = 2, n_neg = 3, n_zer = 1'
  IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
  CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
  CALL find_dyn_cell1(package(pack_index)%pos, cur_pgi)
  CALL change_cell(pack_index, cur_pgi)
 ! ELSE IF(n_par > 1) THEN
 !  
 ELSE IF(n_pos == 3 .and. n_neg == 3) THEN
  change_of_cell = .FALSE.
  IF(procout) write(*,*) 'do_rpackage: n_pos == 3, n_neg == 3'
  CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
 ELSE IF(n_pos == 4 .and. n_neg == 2) THEN
  IF(procout) write(*,*) 'do_rpackage: n_pos == 4 and n_neg == 2'
  cur_pos = package(pack_index)%pos
  cur_pgi = package(pack_index)%cell_numb
  CALL find_dyn_cell1(cur_pos, next_cell)
  CALL change_cell(pack_index, next_cell)
  ! write(25,*) cur_pos, dyn_cell(next_cell)%corner, dyn_cell(next_cell)%upcorner
  ! write(26,*) cur_pos, package(pack_index)%dir
  ! STOP 'do_rpackage: testing'
 ELSE
  change_of_cell = .FALSE.
  IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
  CALL move_package(pack_index, e_dist, next_cell, change_of_cell)
 END IF
 ! if the e_dist == 0 then the packets was re-emmited in the photosphere and no line interaction is allowed
 IF(e_dist > 0.D0) THEN
  IF(procout) write(*,*) 'do_rpackage: move_package, e_dist = ', e_dist, ' > 0.D0'
  CALL update_estimators(pack_index, e_dist)
  CALL do_rpackage_event(pack_index, event, actirrates)
  pack_pos = package(pack_index)%pos
  cur_pgi = package(pack_index)%cell_numb
  cur_corner = dyn_cell(cur_pgi)%corner
  cur_width = dyn_cell(cur_pgi)%upcorner - dyn_cell(cur_pgi)%corner
  ! correction of position
  ! IF((pack_pos(ind_x) < cur_corner(ind_x) + epsilon0 .and.  pack_pos(ind_x) > cur_corner(ind_x) - epsilon0) &
  !  .or. (pack_pos(ind_y) < cur_corner(ind_y) + epsilon0 .and.  pack_pos(ind_y) > cur_corner(ind_y) - epsilon0) .or. &
  !   (pack_pos(ind_z) < cur_corner(ind_z) + epsilon0 .and.  pack_pos(ind_z) > cur_corner(ind_z) - epsilon0) .or.&
  !   (pack_pos(ind_x) < cur_corner(ind_x) + width(ind_x) + epsilon0 &
  !    .and.  pack_pos(ind_x) > cur_corner(ind_x) + width(ind_x) - epsilon0) .or.&
  !   (pack_pos(ind_y) < cur_corner(ind_y) + width(ind_y) + epsilon0 &
  !    .and.  pack_pos(ind_y) > cur_corner(ind_y) + width(ind_y) - epsilon0) .or.&
  !   (pack_pos(ind_z) < cur_corner(ind_z) + width(ind_z) + epsilon0 &
  !    .and.  pack_pos(ind_z) > cur_corner(ind_z) + width(ind_z) - epsilon0)) THEN
  !   IF(procout) write(*,*) 'do_rpackage: correction of position, moving for cur_dist = ', cur_dist
  !   change_of_cell = .FALSE.
  !   cur_dist = 1.D2
  !   CALL move_package(pack_index, cur_dist, next_cell, change_of_cell)
  !   cur_pos = package(pack_index)%pos
  !   CALL find_dyn_cell1(cur_pos, cur_pgi)
  !   CALL change_cell(pack_index, cur_pgi)
  ! END IF
  IF(procout) write(*,*) 'do_rpackage: move_package, e_dist = 0.D0'
  CALL bound_dist(pack_index, cur_pgi, cell_dist, n_pos, n_neg, n_zer, n_par)
 END IF
!_______________________________________________________________________________________________
! event_dist > cell_dist
ELSE IF(e_dist > cell_dist) THEN
 IF(procout) write(*,*) 'do_rpackage: e_dist > cell_dist'
 IF(cell_dist > 0.D0) THEN
  IF(procout) write(*,*) 'do_rpackage: cell_dist > 0'
  IF((n_pos == 3 .and. n_neg == 3) .or. (n_pos == 3 .and. n_neg == 2 .and. n_zer == 1)) THEN 
   ! Move package from the curent position for the cell_dist
   change_of_cell = .TRUE.
   IF(procout) write(*,*) 'do_rpackage: cell_dist > epsilon0 or (cell_dist > 0 and n_pos == 3 and',&
   ' n_neg == 3)'
   IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
   CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
   CALL update_estimators(pack_index, cell_dist)
  !_______________________________________________________________________________________________
  ! ELSE IF(n_pos == 2 .and. n_neg == 3 .and. n_zer == 1) THEN 
  !_______________________________________________________________________________________________
  ELSE IF(cell_dist < epsilon0 .and. n_pos == 4 .and. n_neg == 2) THEN 
   IF(procout) write(*,*) 'do_rpackage: e_dist > cell_dist, cell_dist > 0, cell_dist < epsilon0 ', &
    ' n_pos == 4, n_neg == 2'
   change_of_cell = .FALSE.
   CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
  ELSE IF(n_pos == 2 .and. n_neg == 3 .and. n_zer == 1) THEN 
   IF(procout) write(*,*) 'do_rpackage: n_pos == 2 and n_neg == 3'
   cur_pos = package(pack_index)%pos
   CALL find_dyn_cell1(cur_pos, next_cell)
   cur_corner = dyn_cell(cur_pgi)%corner
   ! CALL change_cell(pack_index, next_cell)
   IF(cur_pos(ind_x) == cur_corner(ind_x)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(negx)
   ELSE IF(cur_pos(ind_y) == cur_corner(ind_y)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(negy)
   ELSE IF(cur_pos(ind_z) == cur_corner(ind_z)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(negz)
   END IF
   CALL change_cell(pack_index, next_cell)

   ! STOP 'do_rpackage: testing'
  ELSE IF(n_pos == 2 .and. n_neg == 2 .and. n_par == 2) THEN
   IF(procout) write(*,*) 'do_rpackage: n_pos = 2 and n_neg = 2 and n_par = 2'
   ! write(50,*) package(pack_index)%pos/R_star, cur_corner/R_star, cur_width/R_star, package(pack_index)%dir
   ! STOP 'do_rpackage: testing'
   change_of_cell = .TRUE.
   CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
  ELSE
   change_of_cell = .TRUE.
   IF(procout) write(*,*) 'do_rpackage: cell_dist > bound_dist and n_par == 2'
   IF(procout) write(*,*) 'do_rpackage: move_package, change_of_cell = ', change_of_cell
   CALL move_package(pack_index, cell_dist, next_cell, change_of_cell)
   CALL update_estimators(pack_index, cell_dist)
   ! write(*,*) 'do_rpackage: n_pos = ', n_pos, ' n_neg = ', n_neg, ' n_zer = ', n_zer, ' n_par = ', n_par
   ! STOP 'do_rpackage: testing'
  END IF
 ELSE IF(cell_dist == 0.D0) THEN
  IF(procout) write(*,*) 'do_rpackage: cell_dist == 0'
 !_______________________________________________________________________________________________
  IF(n_pos == 1 .and. n_neg == 4 .and. n_zer == 1) THEN
   IF(procout) write(*,*) 'do_rpackage: cell_dist == 0 and n_pos == 1 and n_neg == 4 and n_zer == 1'
   IF(procout) write(*,*) 'do_rpackage: setting up new propGrid cell number to: ', cur_pgi
   cur_pos = package(pack_index)%pos
   cur_corner = dyn_cell(cur_pgi)%corner
   cur_upcorner = dyn_cell(cur_pgi)%upcorner
   cur_width = dyn_cell(cur_pgi)%upcorner - dyn_cell(cur_pgi)%corner
   CALL find_dyn_cell1(cur_pos, cur_pgi)
   write(*,*) 'do_rpackage cell starting = ', dyn_cell(cur_pgi)%corner
   write(*,*) 'do_rpackage packet pos = ', package(pack_index)%pos
   write(*,*) 'do_rpackage cell ending = ', dyn_cell(cur_pgi)%upcorner
   IF(cur_pos(ind_x) == cur_corner(ind_x)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(negx)
   ELSE IF(cur_pos(ind_y) == cur_corner(ind_y)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(negy)
   ELSE IF(cur_pos(ind_z) == cur_corner(ind_z)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(negz)
   ELSE IF(cur_pos(ind_x) == cur_upcorner(ind_x)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(posx)
   ELSE IF(cur_pos(ind_y) == cur_upcorner(ind_y)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(posy)
   ELSE IF(cur_pos(ind_z) == cur_upcorner(ind_z)) THEN
    next_cell = dyn_cell(cur_pgi)%neighbor(posz)
   END IF
   CALL change_cell(pack_index, next_cell)
  END IF
 ELSE IF(cell_dist < 0.D0) THEN
  IF(procout) write(*,*) 'do_rpackage: cell_dist < 0'
  IF(((next_cross <= 6) .and. (next_cross >=1))) THEN
   ! next_cell = dyn_cell(cur_pgi)%neighbor(package(pack_index)%next_cross)
   IF(procout) write(*,*) 'do_rpackage: cell_dist = ', cell_dist, ' next_cross = ', next_cross
   IF(procout) write(*,*) 'do_rpackage: only changing the cell'
   next_cross = package(pack_index)%next_cross
   cur_corner = dyn_cell(next_cell)%corner
   cur_upcorner = dyn_cell(next_cell)%upcorner
   cur_width = dyn_cell(next_cell)%width
   IF(next_cross == posx) THEN
    package(pack_index)%pos(ind_x) = cur_upcorner(ind_x)
   ELSE IF(next_cross == negx) THEN
    package(pack_index)%pos(ind_x) = cur_corner(ind_x)
   ELSE IF(next_cross == posy) THEN
    package(pack_index)%pos(ind_y) = cur_upcorner(ind_y)
   ELSE IF(next_cross == negy) THEN
    package(pack_index)%pos(ind_y) = cur_corner(ind_y)
   ELSE IF(next_cross == posz) THEN
    package(pack_index)%pos(ind_z) = cur_upcorner(ind_z)
   ELSE IF(next_cross == negz) THEN
    package(pack_index)%pos(ind_z) = cur_corner(ind_z)
   END IF
   cur_corner = dyn_cell(next_cell)%corner
   cur_width = dyn_cell(next_cell)%width
   CALL find_dyn_cell1(package(pack_index)%pos, cur_pgi)
   CALL change_cell(pack_index, cur_pgi)
  !_______________________________________________________________________________________________
  ELSE IF(next_cross > 6) THEN
   IF(procout) write(*,*) 'do_rpackage: cell_dist < 0 and next_cross > 6'
   cur_corner = dyn_cell(cur_pgi)%corner
   cur_pos = package(pack_index)%pos
   delta_r = (cur_corner - cur_pos)/NORM2(cur_corner - cur_pos)
   ! write(*,*) sign(1.D0, package(pack_index)%dir(ind_x))
   ! write(*,*) sign(1.D0, package(pack_index)%dir(ind_y))
   ! write(*,*) sign(1.D0, package(pack_index)%dir(ind_z))
   if(delta_r(ind_x) == 0.D0) delta_r(ind_x) = 1.D10
   if(delta_r(ind_y) == 0.D0) delta_r(ind_y) = 1.D10
   if(delta_r(ind_z) == 0.D0) delta_r(ind_z) = 1.D10
   IF(next_cross == edyz) THEN                  
    ! move packet to with the vector (0,1,1)
    package(pack_index)%pos = package(pack_index)%pos + (/0.D0,delta_r(ind_y),delta_r(ind_z)/)
    IF(procout) write(*,*) 'do_rpackage: move packet', pack_index, ' with the vector delta_r = ', delta_r
   ELSE IF(next_cross == edxz) THEN
    ! move packet to with the vector (1,0,1)    
    package(pack_index)%pos = package(pack_index)%pos + (/delta_r(ind_x),0.D0,delta_r(ind_z)/)
    IF(procout) write(*,*) 'do_rpackage: move packet', pack_index, ' with the vector delta_r = ', delta_r
   ELSE IF(next_cross == edxy) THEN             
    ! move packet to with the vector (0,1,1)
    package(pack_index)%pos = package(pack_index)%pos + (/delta_r(ind_x),delta_r(ind_y),0.D0/)
    IF(procout) write(*,*) 'do_rpackage: move packet', pack_index, ' with the vector delta_r = ', delta_r
   END IF
   CALL find_dyn_cell1(cur_pos, cur_pgi)
   package(pack_index)%cell_numb = cur_pgi
   CALL boundary3(pack_index, cell_dist, next_cell, n_pos, n_neg, n_zer, n_par)
   write(*,*) 'do_rpackage corr: _________________________________________________'
   write(*,*) 'do_rpackage corr: cell starting = ', dyn_cell(pack_index)%corner
   write(*,*) 'do_rpackage corr: packet pos = ', package(pack_index)%pos
   write(*,*) 'do_rpackage corr: cell ending = ', dyn_cell(pack_index)%upcorner
   write(*,*) 'do_rpackage corr: _________________________________________________'
  END IF
 END IF
END IF ! e_dist > b_dist

!_______________________________________________________________________________________________
cur_pgi = package(pack_index)%cell_numb
pos = package(pack_index)%pos
corner = dyn_cell(cur_pgi)%corner
upcorner = dyn_cell(cur_pgi)%upcorner
! DO ind_I = 1, const_dimofspace
!  IF ( (pos(ind_I) < corner(ind_I) .OR. pos(ind_I) > upcorner(ind_I)) .AND. pack_index <= SIZE(package) ) THEN
!  ! IF((pos(ind_I) < corner(ind_I) .OR. pos(ind_I) > upcorner(ind_I)) .and. pack_index <= SIZE(package)) THEN
!   if(pack_index == 21) write(98,*) pos/const_Rsun, dyn_cell(cur_pgi)%corner/const_Rsun, dyn_cell(cur_pgi)%width/const_Rsun
!   CALL correction_propagation(ind_I, pack_index, cur_pgi)
!   if(pack_index == 21) write(98,*) pos/const_Rsun, dyn_cell(cur_pgi)%corner/const_Rsun, dyn_cell(cur_pgi)%width/const_Rsun
!   cur_pgi = package(pack_index)%cell_numb
!   pos = package(pack_index)%pos
!   corner = dyn_cell(cur_pgi)%corner
!   upcorner = dyn_cell(cur_pgi)%upcorner
!  END IF
! END DO

! check if we are in a correct cell
IF(debug == 2) THEN
 cur_pgi = package(pack_index)%cell_numb
 pos = package(pack_index)%pos
 corner = dyn_cell(cur_pgi)%corner
 upcorner = dyn_cell(cur_pgi)%upcorner
 write(*,*) 'do_rpackage II: cell starting = ', corner
 write(*,*) 'do_rpackage II: packet pos = ', package(pack_index)%pos
 write(*,*) 'do_rpackage II: cell ending = ', upcorner
 CALL find_dyn_cell1(pos, pomocna_bunka)
 write(*,*) 'do_rpackage II: pack_index = ', pack_index, ' bunka = ', pomocna_bunka
 write(*,*) 'do_rpackage II: cur_pgi = ', cur_pgi
 write(*,*) 'do_rpackage II: direction = ', package(pack_index)%dir
 write(*,*) 'do_rpackage II: active = ', package(pack_index)%active
 DO ind_I = 1, const_dimofspace
  IF(((pos(ind_I) < corner(ind_I) - mininum) .OR. (pos(ind_I) > upcorner(ind_I) + mininum)) &
   .and. pack_index <= SIZE(package)) THEN
   write(*,*) 'do_rpackage i: r - min = ', (corner(ind_I) - mininum), ' r + w + min = ', &
    (upcorner(ind_I) + mininum)
   ! CALL find_dyn_cell1(pos, pomocna_bunka)
   ! write(*,*) 'do_rpackage: pos - cell = ', 
   ! write(*,*) 'do_rpackage: skutecna bunka = ', pomocna_bunka
   write(*,*) 'do_rpackage II: _________________________________________________'
   write(*,*) 'do_rpackage II: cell starting = ', dyn_cell(pomocna_bunka)%corner
   write(*,*) 'do_rpackage II: packet pos = ', package(pack_index)%pos
   write(*,*) 'do_rpackage II: cell ending = ', dyn_cell(pomocna_bunka)%upcorner
   write(*,*) 'do_rpackage II: _________________________________________________'
   write(4,*) pos, dyn_cell(cur_pgi)%corner, dyn_cell(cur_pgi)%width, get_package_model_index(pack_index)
   write(4,*) pos, corner, upcorner, get_package_model_index(pack_index)
   write(*,*) 'I = ', ind_I
   ! write(50,*) package(pack_index)%pos/R_star, corner/R_star, width/R_star, package(pack_index)%dir
   STOP 'do_rpackage: packet is not located inside the propagation cell'
  END IF
 END DO
END IF

  
END SUBROUTINE do_rpackage
