! this sbr calculates a dynamics of v-packets
! v-packets can exist outside the propGrid, however it muse be calculated whether this packet flyes into the computational domain
SUBROUTINE do_vpackage(cur_vpackage)


USE types
USE constants

IMPLICIT NONE

INTEGER                                         :: cur_vpackage

DOUBLE PRECISION, DIMENSION(3)                  :: pos, dir
DOUBLE PRECISION, DIMENSION(6)                  :: tecka, poser
DOUBLE PRECISION, DIMENSION(3)                  :: dist
DOUBLE PRECISION, DIMENSION(3)                  :: cross_point
DOUBLE PRECISION, DIMENSION(3)                  :: cur_pos, cur_dir
INTEGER                                         :: cur_poser

DOUBLE PRECISION                                :: A, dummy, dummy2
DOUBLE PRECISION                                :: cur_tecka
DOUBLE PRECISION, PARAMETER                     :: largeNumber = 1.D90
INTEGER                                         :: I, J
INTEGER                                         :: next_cross

DOUBLE PRECISION                                :: D, L_star
INTEGER                                         :: ind_cell_numb, virt_pack_index

LOGICAL                                         :: seeking

L_star = 4.D0*pi*(R_star)**2*sigma*T_eff**4

pos = vpackage(cur_vpackage)%pos
dir = vpackage(cur_vpackage)%dir

! poser: an array of directions
poser = (/ posx, negx, posy, negy, posz, negz /)

! tecka: an array of distances
! package is outside the propGrid
IF(pos(1) < xmin .or. pos(1) > xmax .or. pos(2) < ymin .or. pos(2) > ymax .or. &
 pos(3) < zmin .or. pos(3) > zmax) THEN
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

 write(*,*) 'do_vpackage: surova, tecka = ', tecka
 
 dist(:) = 1.D99
 next_cross = -99

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
 write(*,*) 'do_vpackage: ordered, tecka = ', tecka
 write(*,*) 'do_vpackage: ordered, poser = ', poser

 ! now going through value by value, we are seeking for the cross point with the propGrid
 cur_tecka = 0
 DO WHILE(seeking)
  cur_tecka = cur_tecka + 1
  cur_poser = poser(cur_tecka)
  cur_pos = pos + cur_poser * dir
  
  IF(cur_poser == posx .or. cur_poser == negx) THEN
   IF(cur_pos(2) > ymin .and. cur_pos(2) < ymax .and. &
    cur_pos(3) > zmin .and. cur_pos(3) < zmax) THEN
     vpackage(cur_vpackage)%pos = cur_pos
     seeking = .false.
   END IF
  END IF
  IF(cur_poser == posy .or. cur_poser == negy) THEN
   IF(cur_pos(1) > xmin .and. cur_pos(1) < xmax .and. &
    cur_pos(3) > zmin .and. cur_pos(3) < zmax) THEN
     vpackage(cur_vpackage)%pos = cur_pos
     seeking = .false.
   END IF
  END IF
  IF(cur_poser == posz .or. cur_poser == negz) THEN
   IF(cur_pos(3) > ymin .and. cur_pos(3) < ymax .and. &
    cur_pos(1) > xmin .and. cur_pos(1) < xmax) THEN
     vpackage(cur_vpackage)%pos = cur_pos
     seeking = .false.
   END IF
  END IF

  ! nothing was found, the packet flies totally out of the propGrid
  IF(cur_tecka > 6) THEN
   seeking = .false.
   vpackage(cur_vpackage)%active = -99
  END IF

 END DO
  write(*,*) 'do_vpackage: pos = ', vpackage(cur_vpackage)%pos/R_inf
 
! the packet is inside the propagation grid
END IF

IF(pos(1) > xmin .or. pos(1) < xmax .or. pos(2) > ymin .or. pos(2) < ymax .or. &
 pos(3) > zmin .or. pos(3) < zmax) THEN
 ! we generate a packet which will be propagatet through the propGrid
 write(*,*) 'do_vpackage: a generation of a package'
 virt_pack_index = SIZE(package) - n_add_pack + 3
 cur_pos = vpackage(cur_vpackage)%pos
 cur_dir = vpackage(cur_vpackage)%dir
 package(virt_pack_index)%pos = cur_pos
 package(virt_pack_index)%dir = cur_dir
 package(virt_pack_index)%freq_rf = vpackage(cur_vpackage)%freq_rf
 CALL find_dyn_cell1(cur_pos,ind_cell_numb)
 package(virt_pack_index)%cell_numb = ind_cell_numb
 package(virt_pack_index)%active     = 1
 package(virt_pack_index)%typ        = type_rpkt
 package(virt_pack_index)%n_interactions = 0
 package(virt_pack_index)%next_cross = NONE
 package(virt_pack_index)%e_rf = 1.e10  
 IF(norm2(cur_pos) < R_inf) THEN
  CALL doppler_factor(virt_pack_index, cur_pos, cur_dir, D)
 ELSE
  D = 1.D0
 END IF
 package(virt_pack_index)%freq_cmf = package(virt_pack_index)%freq_rf * D 
 package(virt_pack_index)%e_cmf    = package(virt_pack_index)%e_rf * D  
 package(virt_pack_index)%last_line = no_line
 package(virt_pack_index)%delta_s = 0.D0
 DO  WHILE (package(virt_pack_index)%active .EQ. 1)
  IF (package(virt_pack_index)%typ .EQ. type_rpkt) THEN
   ! IF(package(virt_pack_index)%n_interactions .EQ. 1000000) THEN
   !  ! write(*,*) 'package ', virt_pack_index, ' interacted for 2000000 times and will be destroyed...'
   !  package(virt_pack_index)%active = 0
   ! END IF
     ! If the packet is of type rpkt, it represents a photon. So it needs to be propagated.
     CALL do_rpackage(virt_pack_index)
  ELSE IF (package(virt_pack_index)%typ .EQ. type_kpkt) THEN 
     ! If the packet is of type kpkt, it represents thermal kinetic energy.
     ! Sample all possible cooling processes and randomly select one of them
      CALL do_kpackage(virt_pack_index)
  ELSE IF (package(virt_pack_index)%typ .EQ. type_ipkt) THEN 
     ! If the packet is of type ipkt, it represents atomic internal energy (excitation/ionization).
     ! Calculate all transition probabilities and randomly select one of them (macro-atom formalism)
     CALL do_ipackage(virt_pack_index)
  ELSE IF (package(virt_pack_index)%typ == type_dpkt) THEN
   ! this state corresponds to diffusive approximation
   CALL do_dpackage(virt_pack_index)
  ELSE
     STOP 'ERROR unknown package typ'
  END IF
 END DO







END IF









END SUBROUTINE do_vpackage
