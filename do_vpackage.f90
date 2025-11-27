! this sbr calculates a dynamics of v-packets
! v-packets can exist outside the propGrid, however it muse be calculated whether
! this packet flyes into the computational domain
! 
! INPUT: cur_vpackage(INT): index of a package
! OUTPUT: NONE
! 
SUBROUTINE do_vpackage(cur_vpackage)


USE types
USE constants

IMPLICIT NONE

INTEGER                                         :: cur_vpackage

DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: pos, dir
DOUBLE PRECISION, DIMENSION(6)                  :: tecka
DOUBLE PRECISION, DIMENSION(3)                  :: dist
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: cur_pos, cur_dir
INTEGER, DIMENSION(6)                           :: poser
INTEGER                                         :: cur_poser, cur_pgi, dummy2

DOUBLE PRECISION                                :: A, dummy, cur_tecka
INTEGER                                         :: cur_index
DOUBLE PRECISION, PARAMETER                     :: largeNumber = 1.D90
INTEGER                                         :: ind_I, ind_J
INTEGER                                         :: next_cross

DOUBLE PRECISION                                :: doppler_D! , L_star
INTEGER                                         :: ind_cell_numb

LOGICAL                                         :: seeking, procout=.false.
DOUBLE PRECISION, PARAMETER                     :: delta_tecka = 1.D2
DOUBLE PRECISION                                :: crossing_x, crossing_y, crossing_z

! L_star = 4.D0*pi*(R_star)**2*const_stefbolz*T_eff**4


! poser: an array of directions
poser = (/ posx, negx, posy, negy, posz, negz /)

! tecka: an array of distances
! package is outside the propGrid
pos = package(cur_vpackage)%pos
dir = package(cur_vpackage)%dir
IF(pos(ind_x) < xmin .or. pos(ind_x) > xmax .or. pos(ind_y) < ymin .or. pos(ind_y) > ymax .or. &
 pos(ind_z) < zmin .or. pos(ind_z) > zmax) THEN
 IF(procout) write(*,*) 'do_vpackage: the packet is outside the propGrid'
 !________________________________________________________________________
 ! we have to find out, whether the packet crosses the compuational domain
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! there could be a problem with indeces !!!
 ! must be tested!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF(dir(ind_x) /= 0) THEN
  tecka(posx) = (xmax - pos(ind_x))/dir(ind_x)
  tecka(negx) = (xmin - pos(ind_x))/dir(ind_x)
 ELSE
  tecka(posx) = largeNumber
  tecka(negx) = largeNumber
 END IF
 IF(dir(ind_y) /= 0) THEN
  tecka(posy) = (ymax - pos(ind_y))/dir(ind_y)
  tecka(negy) = (ymin - pos(ind_y))/dir(ind_y)
 ELSE
  tecka(posy) = largeNumber
  tecka(negy) = largeNumber
 END IF
 IF(dir(ind_z) /= 0) THEN
  tecka(posz) = (zmax - pos(ind_z))/dir(ind_z)
  tecka(negz) = (zmin - pos(ind_z))/dir(ind_z)
 ELSE
  tecka(posz) = largeNumber
  tecka(negz) = largeNumber
 END IF
 
 dist(:) = 1.D99
 next_cross = -99

 !________________________________________________
 ! order distances from the lowerst to the highest
 DO ind_J = 2, 6
  ind_I = ind_J - 1
  
  A = tecka(ind_J)

  DO WHILE(ind_I >= 1)
   IF(tecka(ind_I) > A) THEN
    dummy = tecka(ind_I + 1)
    dummy2 = poser(ind_I + 1)
    tecka(ind_I + 1) = tecka(ind_I)
    poser(ind_I + 1) = poser(ind_I)
    tecka(ind_I) = dummy
    poser(ind_I) = dummy2
   END IF
    ind_I = ind_I - 1
  END DO
 END DO
 ! write(*,*) 'do_vpackage: ordered, tecka = ', tecka
 ! write(*,*) 'do_vpackage: ordered, poser = ', poser

 !_______________________________________________________________________________________
 ! now going through value by value, we are seeking for the cross point with the propGrid
 cur_index = 0
 seeking = .true.
 DO WHILE(seeking)
  cur_index = cur_index + 1
  cur_tecka = tecka(cur_index)
  cur_poser = poser(cur_index)
  cur_pos = pos + (cur_tecka+delta_tecka) * dir
  CALL find_dyn_cell1(cur_pos, cur_pgi)
  ! write(*,*) 'do_vpackage: cur_pgi = ', cur_pgi
  package(cur_vpackage)%cell_numb = cur_pgi
  ! STOP 'do_vpackage: testing'
  ! write(*,*) 'do_vpackage: dir = ', dir
  ! write(*,*) 'do_vpackage: cur_index = ', cur_index, ' cur_tecka = ', cur_tecka
  ! write(*,*) 'do_vpackage: cur_poser = ', cur_poser
  ! write(*,*) 'do_vpackage: position = ', cur_pos(ind_x)/xmin, cur_pos(ind_y)/ymin, cur_pos(ind_z)/zmin
  
  IF(cur_poser == posx .or. cur_poser == negx) THEN
   IF(cur_pos(ind_y) > ymin .and. cur_pos(ind_y) < ymax .and. &
    cur_pos(ind_z) > zmin .and. cur_pos(ind_z) < zmax) THEN
     package(cur_vpackage)%pos = cur_pos
     seeking = .false.
    IF(procout) write(*,*) 'do_vpackage: posx/negx was chosen'
    ! write(64,*) dir, cur_pos(:)/xmax
   END IF
  END IF
  IF(cur_poser == posy .or. cur_poser == negy) THEN
   IF(cur_pos(ind_x) > xmin .and. cur_pos(ind_x) < xmax .and. &
    cur_pos(ind_z) > zmin .and. cur_pos(ind_z) < zmax) THEN
     ! write(*,*) 'do_vpackage: setting a new position for the vpackage'
     package(cur_vpackage)%pos = cur_pos
     seeking = .false.
    IF(procout) write(*,*) 'do_vpackage: posy/negy was chosen'
    ! write(64,*) dir, cur_pos(:)/xmax
   END IF
  END IF
  IF(cur_poser == posz .or. cur_poser == negz) THEN
   IF(cur_pos(ind_z) > ymin .and. cur_pos(ind_z) < ymax .and. &
    cur_pos(ind_x) > xmin .and. cur_pos(ind_x) < xmax) THEN
     package(cur_vpackage)%pos = cur_pos
     seeking = .false.
    IF(procout) write(*,*) 'do_vpackage: posz/negz was chosen'
    ! write(64,*) dir, cur_pos(:)/xmax
   END IF
  END IF

  ! write(*,*) 'do_vpackage: seeking = ', seeking

  ! nothing was found, the packet flies totally out of the propGrid
  IF(cur_index == 6 .and. seeking .EQV. .true.) THEN
   seeking = .false.
   package(cur_vpackage)%active = -99
   package(cur_vpackage)%cell_numb = cur_pgi
   IF(procout) write(*,*) 'do_vpackage: nothing was found'
  END IF

 END DO ! while seeking
 ! write(*,*) 'do_vpackage: active = ', package(cur_vpackage)%active
 ! write(*,*) 'do_vpackage: pos = ', package(cur_vpackage)%pos(ind_y)/ymin
 
END IF
 pos = package(cur_vpackage)%pos
!_______________________________________________________________________________________
!______________________ INSIDE THE PROPGRID ____________________________________________
!_______________________________________________________________________________________
! the packet is inside the propagation grid
IF(debug == 2) THEN
 write(*,*) 'do_vpackage: pos = ', pos, ' y/ymin = ',pos(ind_y)/ymin, 'y/ymax = ', pos(ind_y)/ymax
END IF

IF(pos(ind_x) >= xmin .and. pos(ind_x) <= xmax .and. pos(ind_y) >= ymin .and. pos(ind_y) <= ymax .and. &
  pos(ind_z) >= zmin .and. pos(ind_z) <= zmax) THEN
 ! we generate a packet which will be propagatet through the propGrid
 ! transcription of the virtual packet properties to the standard packet properties
 IF(procout) write(*,*) 'do_vpackage: the packet is inside the propGrid'
 IF(debug == 2) THEN
  write(*,*) 'do_vpackage: pos/R_inf = ', pos/R_inf
 END IF
 cur_pos = package(cur_vpackage)%pos

 IF(debug == 2) THEN
  write(*,*) 'do_vpackage: init pos = ', package(cur_vpackage)%pos/xmax
 END IF
 if(procout) write(*,*) 'do_vpackage: init pos = ', package(cur_vpackage)%pos/xmax

 cur_dir = package(cur_vpackage)%dir
 CALL find_dyn_cell1(cur_pos,ind_cell_numb)
 package(cur_vpackage)%cell_numb = ind_cell_numb

 IF(ind_cell_numb <= 0) THEN
  write(*,*) 'do_vpackage: the propGrid index cell number = ', ind_cell_numb, ' is lower than zero!'
  STOP
 END IF

 IF(norm2(cur_pos) < R_inf) THEN
  ! CALL doppler_factor(cur_vpackage, cur_pos, cur_dir, doppler_D)
  CALL doppler_factor(cur_vpackage, doppler_D)
 ELSE
  doppler_D = 1.D0
 END IF

 package(cur_vpackage)%freq_cmf = package(cur_vpackage)%freq_rf * doppler_D 
 package(cur_vpackage)%e_cmf    = package(cur_vpackage)%e_rf * doppler_D


 ! propagation of the packet
 ! write(*,*) 'do_vpackage: cur_vpackage = ', cur_vpackage
 CALL packet_dynamics(cur_vpackage)


 IF(debug == 2) THEN
  write(*,*) 'do_vpackage: active = ', package(cur_vpackage)%active
  write(*,*) 'do_vpackage: end pos = ', package(cur_vpackage)%pos/xmax
 END IF

END IF ! packet is inside the propGrid




END SUBROUTINE do_vpackage
