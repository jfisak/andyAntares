! this sbr calculates a dynamics of v-packets
! v-packets can exist outside the propGrid, however it muse be calculated whether this packet flyes into the computational domain
SUBROUTINE do_vpackage(cur_vpackage)


USE types
USE constants

IMPLICIT NONE

INTEGER                                         :: cur_vpackage

DOUBLE PRECISION, DIMENSION(3)                  :: pos, dir
DOUBLE PRECISION, DIMENSION(6)                  :: tecka
DOUBLE PRECISION, DIMENSION(3)                  :: dist
DOUBLE PRECISION, DIMENSION(3)                  :: cur_pos, cur_dir
INTEGER, DIMENSION(6)                           :: poser
INTEGER                                         :: cur_poser, cur_pgi, dummy2

DOUBLE PRECISION                                :: A, dummy, cur_tecka
INTEGER                                         :: cur_index
DOUBLE PRECISION, PARAMETER                     :: largeNumber = 1.D90
INTEGER                                         :: I, J
INTEGER                                         :: next_cross

DOUBLE PRECISION                                :: D! , L_star
INTEGER                                         :: ind_cell_numb

LOGICAL                                         :: seeking
DOUBLE PRECISION, PARAMETER                     :: delta_tecka = 1.D2

! L_star = 4.D0*pi*(R_star)**2*sigma*T_eff**4


! poser: an array of directions
poser = (/ posx, negx, posy, negy, posz, negz /)

! tecka: an array of distances
! package is outside the propGrid
pos = package(cur_vpackage)%pos
dir = package(cur_vpackage)%dir
IF(pos(1) < xmin .or. pos(1) > xmax .or. pos(2) < ymin .or. pos(2) > ymax .or. &
 pos(3) < zmin .or. pos(3) > zmax) THEN
 !________________________________________________________________________
 ! we have to find out, whether the packet crosses the compuational domain
 IF(dir(1) /= 0) THEN
  tecka(1) = (xmin - pos(1))/dir(1)
  tecka(2) = (xmax - pos(1))/dir(1)
 ELSE
  tecka(1) = largeNumber
  tecka(2) = largeNumber
 END IF
 IF(dir(2) /= 0) THEN
  tecka(3) = (ymin - pos(2))/dir(2)
  tecka(4) = (ymax - pos(2))/dir(2)
 ELSE
  tecka(3) = largeNumber
  tecka(4) = largeNumber
 END IF
 IF(dir(3) /= 0) THEN
  tecka(5) = (zmin - pos(3))/dir(3)
  tecka(6) = (zmax - pos(3))/dir(3)
 ELSE
  tecka(5) = largeNumber
  tecka(6) = largeNumber
 END IF
 
 dist(:) = 1.D99
 next_cross = -99

 !________________________________________________
 ! order distances from the lowerst to the highest
 DO J = 2, 6
  I = J - 1
  
  A = tecka(J)

  DO WHILE(I >= 1)
   IF(tecka(I) > A) THEN
    dummy = tecka(I + 1)
    dummy2 = poser(I + 1)
    tecka(I + 1) = tecka(I)
    poser(I + 1) = poser(I)
    tecka(I) = dummy
    poser(I) = dummy2
   END IF
    I = I - 1
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
  ! write(*,*) 'do_vpackage: position = ', cur_pos(1)/xmin, cur_pos(2)/ymin, cur_pos(3)/zmin
  
  IF(cur_poser == posx .or. cur_poser == negx) THEN
   IF(cur_pos(2) > ymin .and. cur_pos(2) < ymax .and. &
    cur_pos(3) > zmin .and. cur_pos(3) < zmax) THEN
     package(cur_vpackage)%pos = cur_pos
     seeking = .false.
   END IF
  END IF
  IF(cur_poser == posy .or. cur_poser == negy) THEN
   IF(cur_pos(1) > xmin .and. cur_pos(1) < xmax .and. &
    cur_pos(3) > zmin .and. cur_pos(3) < zmax) THEN
     write(*,*) 'do_vpackage: setting a new position for the vpackage'
     package(cur_vpackage)%pos = cur_pos
     seeking = .false.
   END IF
  END IF
  IF(cur_poser == posz .or. cur_poser == negz) THEN
   IF(cur_pos(3) > ymin .and. cur_pos(3) < ymax .and. &
    cur_pos(1) > xmin .and. cur_pos(1) < xmax) THEN
     package(cur_vpackage)%pos = cur_pos
     seeking = .false.
   END IF
  END IF

  ! write(*,*) 'do_vpackage: seeking = ', seeking

  ! nothing was found, the packet flies totally out of the propGrid
  IF(cur_index == 6 .and. seeking .EQV. .true.) THEN
   seeking = .false.
   package(cur_vpackage)%active = -99
   package(cur_vpackage)%cell_numb = cur_pgi
  END IF

 END DO ! while seeking
 ! write(*,*) 'do_vpackage: active = ', package(cur_vpackage)%active
 ! write(*,*) 'do_vpackage: pos = ', package(cur_vpackage)%pos(2)/ymin
 
END IF
 pos = package(cur_vpackage)%pos
!_______________________________________________________________________________________
!______________________ INSIDE THE PROPGRID ____________________________________________
!_______________________________________________________________________________________
! the packet is inside the propagation grid
IF(debug == 2) THEN
 write(*,*) 'do_vpackage: pos = ', pos, ' y/ymin = ',pos(2)/ymin, 'y/ymax = ', pos(2)/ymax
END IF

IF(pos(1) >= xmin .and. pos(1) <= xmax .and. pos(2) >= ymin .and. pos(2) <= ymax .and. &
  pos(3) >= zmin .and. pos(3) <= zmax) THEN
 ! we generate a packet which will be propagatet through the propGrid
 ! transcription of the virtual packet properties to the standard packet properties
 IF(debug == 2) THEN
  write(*,*) 'do_vpackage: pos/R_inf = ', pos/R_inf
 END IF
 cur_pos = package(cur_vpackage)%pos

 IF(debug == 2) THEN
  write(*,*) 'do_vpackage: init pos = ', package(cur_vpackage)%pos/xmax
 END IF

 cur_dir = package(cur_vpackage)%dir
 CALL find_dyn_cell1(cur_pos,ind_cell_numb)
 package(cur_vpackage)%cell_numb = ind_cell_numb

 IF(ind_cell_numb <= 0) THEN
  write(*,*) 'do_vpackage: the propGrid index cell number = ', ind_cell_numb, ' is lower than zero!'
  STOP
 END IF

 IF(norm2(cur_pos) < R_inf) THEN
  ! CALL doppler_factor(cur_vpackage, cur_pos, cur_dir, D)
  CALL doppler_factor(cur_vpackage, D)
 ELSE
  D = 1.D0
 END IF

 package(cur_vpackage)%freq_cmf = package(cur_vpackage)%freq_rf * D 
 package(cur_vpackage)%e_cmf    = package(cur_vpackage)%e_rf * D  


 ! propagation of the packet
 ! write(*,*) 'do_vpackage: cur_vpackage = ', cur_vpackage
 CALL packet_dynamics(cur_vpackage)


 IF(debug == 2) THEN
  write(*,*) 'do_vpackage: active = ', package(cur_vpackage)%active
  write(*,*) 'do_vpackage: end pos = ', package(cur_vpackage)%pos/xmax
 END IF

END IF ! packet is inside the propGrid




END SUBROUTINE do_vpackage
