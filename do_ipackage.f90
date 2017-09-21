! this subroutine computes an i-packet dynamics
SUBROUTINE do_ipackage(pack_index)

USE TYPES
USE rates
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index
INTEGER                         :: actual_state
INTEGER                         :: last_line, last_level, last_ion
! is a macro atom active?
INTEGER                         :: active
! number of line transitions
INTEGER                         :: nlns, nluns
INTEGER                         :: element_index, ion_index
INTEGER                         :: I, J, K, line
INTEGER, ALLOCATABLE            :: linetransitions(:), transitions(:), lineuptransitions(:), &
                                   lineradtransitions(:) ! an array of transitions down from last_line
INTEGER                         :: nlevslion
DOUBLE PRECISION                :: actVal
DOUBLE PRECISION                :: ran2, rand
! sum function
DOUBLE PRECISION                :: Z, Ztotal, Zintdownrad, Zraddeexc, Zintuprad, Zrad, Zcoll, &
                                   Zintupcoll, Zintdowncoll, Zdown, Zup, Zintdown, Zintup, &
                                   Zphotrecom, Zintrecombination, Zcollrecom
! B-F internal processes
DOUBLE PRECISION                :: Zphotionup, Zphotiondown, Zcollionup, Zcolliondown
DOUBLE PRECISION                :: Zionization, Zrecombination
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7
DOUBLE PRECISION                :: summ, stat_weight, exci_energy, exci_energy_l, exci_energy_u
! populations
DOUBLE PRECISION                :: act_pop
INTEGER                         :: get_package_model_index, current_mgi
! new frequency
DOUBLE PRECISION                :: new_freq
! Doppler factor
DOUBLE PRECISION                :: D

! define the needed variables
! it is necessary to remember the initial conditions of a macro-atom
last_line = package(pack_index)%last_line
last_ion = linelist(last_line)%indexi
last_level = linelist(last_line)%upper
linelist(last_line)%n_exc = linelist(last_line)%n_exc + 1
element_index = linelist(last_line)%indexe
ion_index = linelist(last_line)%indexi
current_mgi = get_package_model_index(pack_index)
!print*, 'do_ipackage: element_index = ', element_index, 'ion_index = ', ion_index, ' level_index = ', last_level

active = 1
! this is an initial state of the macro-atom
actual_state = last_level
! we will run this loop until the macro atom is deactivated
DO WHILE (active == 1)
 !print*, 'do_ipackage: actual_state = ', actual_state, ' ion_index = ', ion_index
 ! we have to find all possible downward upward transitions
 ! firstly we calculate number of these possible transitions
 ! number of transitions to a lower level
 nlns = 0
 ! number of transitions to a upper level
 nluns = 0
 DO I = 1, ntransitions
  ! we are interested only in the transitions for the given atom
  IF(linelist(I)%indexe == element_index .AND. linelist(I)%indexi == ion_index) THEN
   ! transitions to a lower level
   IF(linelist(I)%upper == actual_state) THEN
    nlns = nlns + 1
   END IF
   ! transitions to a upper level
   IF(linelist(I)%lower == actual_state) THEN
    nluns = nluns + 1
   END IF
  END IF
 END DO
 
 !print*, 'do_ipackage: number of found transitions: ', nlns
 ALLOCATE(linetransitions(nlns), lineuptransitions(nluns))
 ! we will save these possible transitions into an array
 J = 0
 K = 0
 DO I = 1, ntransitions
  IF(linelist(I)%indexe == element_index .AND. linelist(I)%indexi == ion_index) THEN
   ! transitions to a lower level
   IF(linelist(I)%upper == actual_state) THEN
    J = J + 1
    linetransitions(J) = I
   END IF
   ! transitions to a upper level
   IF(linelist(I)%lower == actual_state) THEN
    K = K + 1
    lineuptransitions(K) = I
   END IF
  END IF
 END DO
 
 ! total rates of procedure
 ! 0.) internal downward jump within the current ion
 ALLOCATE(Lma_int_dorad(nlns))
 ! 1.) radiative deexcitation
 ! this allocates only in the first loop, because this field will only remember transitions
 ! from the last_line's upper level
 ! 3.) collisional deexcitation
 ALLOCATE(Lma_int_docoll(nlns))
 ! 2.) internal upward jump within the current ion
 ALLOCATE(Lma_int_uprad(nluns))
 ALLOCATE(Lma_int_upcoll(nluns))
 ALLOCATE(Lma_int_do(nlns), Lma_int_up(nluns))
 Zintdown = 0.D0
 Zintup = 0.D0
 ! allocation of the field for recombination processes
 IF(ion_index > 1) THEN
  nlevslion = SIZE(elements(element_index)%ions(ion_index - 1)%levels)
  !print*, 'element_index = ', element_index, 'ion_index - 1 = ', ion_index - 1, ' nlevslion = ', nlevslion
  ALLOCATE(Lma_recrad(nlevslion), Lma_int_recrad(nlevslion), Lma_reccol(nlevslion), Lma_int_reccol(nlevslion))
 ELSE
  nlevslion = 0
 END IF
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! calculation of the given transition probabilities
 CALL populations(element_index, ion_index, actual_state, current_mgi, act_pop)
 CALL i_radtrans(nlns, linetransitions, nluns, lineuptransitions, act_pop, &
  Zintdownrad, Zintuprad, Zraddeexc)
 CALL i_coltrans(1, pack_index, actual_state, nlns, linetransitions, nluns, lineuptransitions, &
  act_pop, Zintdowncoll, Zintupcoll, Zcoll)
 CALL i_radion(0, element_index, ion_index, actual_state, current_mgi, act_pop, Zphotionup, &
  Zphotiondown, Zphotrecom)
 CALL i_colion(1, element_index, ion_index, actual_state, pack_index, act_pop, Zcollionup, &
  Zcolliondown,Zcollrecom)
 ! total rates of internal donwnward jump
 DO I = 1, nlns
   Lma_int_do(I) = Lma_int_dorad(I) + Lma_int_docoll(I)
 END DO
 Zintdown = Zintdownrad + Zintdowncoll
 ! total rates of internal upward jump
 DO I = 1, nluns
  ! internal jump up
   Lma_int_up(I) = Lma_int_uprad(I) + Lma_int_upcoll(I)
   !print*, 'Lma_int_uprad(I) = ', Lma_int_uprad(I), ' Lma_int_upcoll(I) = ', Lma_int_upcoll(I)
 END DO
 Zintup = Zintuprad + Zintupcoll
! an internal ionization sum
Zionization = Zphotionup + Zcollionup
Zintrecombination = Zphotiondown + Zcolliondown
Zrecombination = Zphotrecom + Zcollrecom
!write(*,*) 'do_ipackage: Zcollrecom = ', Zcollrecom
! the total sum 
Ztotal = Zintdown + Zraddeexc + Zintup + Zcoll + Zionization + Zrecombination + Zintrecombination
! a random number for computation, which process occurs
rand = ran2(idum)
!print*, 'do_ipackage: Ztotal = ', Ztotal
rand = rand * Ztotal
! print*, 'do_ipackage: random number: ', rand, ' Ztotal = ', Ztotal
! these variables are only to the whole line won't be too long
!write(*,*) 'do_ipackage: Zintup = ', Zintup
Z0 = Zintdown
Z1 = Z0 + Zraddeexc
Z2 = Z1 + Zintup
Z3 = Z2 + Zcoll
Z4 = Z3 + Zionization
Z5 = Z4 + Zintrecombination
Z6 = Z5 + Zphotrecom
Z7 = Z6 + Zcollrecom
!print*, 'Zintdown, Zintup, Zraddeexc, Zcoll, Zionization, Zrecombination, Zintrecombination: ', &
!        Zintdown, Zintup, Zraddeexc, Zcoll, Zionization, Zrecombination, Zintrecombination
!write(*,*) 'do_ipackage: Zrecombination = ', Zrecombination, ' Zphotrecom = ', Zphotrecom, ' Zcollrecom = ', Zcollrecom
!write(*,*) 'do_ipackage: Z0 = ', Z0, ' Z1 = ', Z1, ' Z2 = ', Z2, ' Z3 = ', Z3, ' Z4 = ', Z4, ' Z5 = ', Z5, ' Z6 =', Z6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump
! in this case a macro-atom transits into a lower state without an energy emission
IF(rand >= 0.D0 .AND. rand < Z0) THEN
 count_intdownjump = count_intdownjump + 1
! next transition will be an internal downward jump
 summ = 0.D0
 DO I = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lma_int_do(I)) THEN
   actual_state = linelist(linetransitions(I))%lower
   !print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
   EXIT
  END IF
  summ = summ + Lma_int_do(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radiative deexciation
! in this case the packet is now transformed back to the r-packet
! we have to choose a new frequency, which would be calculated randomly from the 
! possible transition last_line -> some lower line
ELSE IF (rand >= Z0 .AND. rand <= Z1) THEN
! next transition will be an radiative deexcitation
 package(pack_index)%typ = type_rpkt
 package(pack_index)%last_line = no_line
 CALL emit_rpackage(pack_index)
 ! now we will calculate new frequency of the packet
 ! we will choose this frequency from the possible radiative transitions
 rand = ran2(idum)
 Zrad = 0.D0
 nlns = 0
 DO K = 1, ntransitions
  ! we are interested only in the transitions for the given atom
  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == linelist(last_line)%indexi) THEN
   ! transitions to a lower level
  ! print*, 'do_ipackage: linelist(K)%upper = ', linelist(K)%upper
   IF(linelist(K)%upper == last_level) THEN
    nlns = nlns + 1
   END IF
  END IF
 END DO
 ALLOCATE(lineradtransitions(nlns), Lma_rad(nlns))
 J = 0
!  ' linelist(last_line)%indexe = ', linelist(last_line)%indexe,&
!  ' linelist(last_line)%indexi = ', linelist(last_line)%indexi, ' last_level = ', last_level, &
!  ' last_line = ', last_line, ' ntransitions = ', ntransitions
 DO K = 1, ntransitions
  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == linelist(last_line)%indexi) THEN
   ! transitions to a lower level
   IF(linelist(K)%upper == last_level) THEN
    J = J + 1
    lineradtransitions(J) = K
    exci_energy_u = &
     elements(element_index)%ions(last_ion)%levels(linelist(K)%upper)%exci_energy
    exci_energy_l = &
     elements(element_index)%ions(last_ion)%levels(linelist(K)%lower)%exci_energy
    stat_weight = &
     elements(element_index)%ions(last_ion)%levels(linelist(K)%lower)%stat_waight
    Lma_rad(J) = linelist(K)%A_ul * stat_weight * (exci_energy_u - exci_energy_l)
    !print*, 'do_ipackage: Lrad(J) = ', Lrad(J)
    Zrad = Zrad + Lma_rad(J)
   END IF
  END IF
 END DO
 rand = rand * Zrad
 summ = 0
 ! looking for the given line
 DO line = 1, nlns
  IF(rand >= summ .AND. rand <= summ + Lma_rad(line)) THEN
   !print*, 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
   !print*, 'raddeexc: summ = ', summ, ' rand = ', rand, ' Zrad = ', Zrad
   ! we found the given cell now we have to compute only a new frequency
   new_freq = linelist(lineradtransitions(line))%freq
   package(pack_index)%freq_cmf = new_freq
   CALL doppler_factor(pack_index, D)
   package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
   package(pack_index)%last_line = no_line
   linelist(lineradtransitions(line))%n_deexc = linelist(lineradtransitions(line))%n_deexc + 1
   IF(linelist(lineradtransitions(line))%lower == linelist(last_line)%lower) THEN
    ! resonant scattering occures
    count_resscattering = count_resscattering + 1
   ELSE
    count_fluorescence = count_fluorescence + 1
   END IF
   EXIT
  END IF
  summ = summ + Lma_rad(line)
 END DO
 active = 0
 DEALLOCATE(Lma_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF (rand >= Z1 .AND. rand <= Z2) THEN
 count_intupjump = count_intupjump + 1
 summ = Z1
 DO I = 1, nluns
  ! we will find the given state
  !write(*,*) 'do_ipackage: summ = ', summ, ' rand = ', rand, ' summ + Lma_int_up = ', summ + Lma_int_up(I)
  IF(rand >= summ .AND. rand < summ + Lma_int_up(I)) THEN
   actual_state = linelist(lineuptransitions(I))%upper
   !print*, 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
   EXIT
  END IF
  summ = summ + Lma_int_up(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional deexcitation
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF(rand >= Z2 .AND. rand <= Z3) THEN
 summ = Z2
 package(pack_index)%last_line = no_line
 !print*, 'pack_index = ', pack_index, ' collisional deexcitation...'
! DO I = 1, nlns
!  ! we will find the given state
!  IF(rand >= summ .AND. rand < summ + Ldowncoll(I)) THEN
!   ! now it will transform into a k-packet, we will have to decide, which k-packet it
!   ! will be
!   CALL col_deexcitation_event(element_index, ion_index, last_line, nlns, linetransitions)
!   EXIT
!  END IF
!  summ = summ + Ldowncoll(I)
! END DO
 !print*, 'collisional deexcitation occures...'
 count_coldeexc = count_coldeexc + 1
 !package(pack_index)%typ = type_kpkt
 ! for now 
 package(pack_index)%typ = type_kpkt
 active = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal photoionization
ELSE IF(rand >= Z3 .AND. rand <= Z4) THEN
 !print*, 'pack_index = ', pack_index, ' internal jump to to the upper ionization state...'
 ion_index = ion_index + 1
 actual_state = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal recombination
ELSE IF(rand >= Z4 .AND. rand <= Z5) THEN
 !print*, 'pack_index = ', pack_index, ' internal jump to to the lower ionization state...'
 ion_index = ion_index - 1
 summ = 0.D0
 rand = ran2(idum)
 rand = rand * Zrecombination
 DO I = 1, nlevslion
  ! we will find the given state
  !print*, 'summ = ', summ, ' summ + L(I) = ', summ + Lma_recrad(I) + Lma_reccol(I)
  IF(rand >= summ .AND. rand < summ + Lma_int_recrad(I) + Lma_int_reccol(I)) THEN
   actual_state = I
!   print*, 'do_ipackage: changing actual state to the state I = ', I
   !print*, 'do_ipackage: packet: ', pack_index, ' internal jump to the lower ionization state...'
   EXIT
  END IF
  summ = summ + Lma_int_recrad(I) + Lma_int_reccol(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radiative recombination
ELSE IF(rand >= Z5 .AND. rand <= Z6) THEN
 !print*, 'do_ipackage: pack_index = ', pack_index, 'radiative recombination'
 package(pack_index)%typ = type_rpkt
 count_rrecombination = count_rrecombination + 1
 active = 0
 ! temporary
 ! the frequency should be sampled from the photion cross section
 count_rrecombination = count_rrecombination + 1
 CALL i_freq_recomb(actual_state, pack_index, act_pop, new_freq)
 package(pack_index)%freq_cmf = new_freq
 CALL doppler_factor(pack_index, D)
 package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional recombination
ELSE IF(rand >= Z6 .AND. rand <= Z7) THEN
 !print*, 'do_ipackage: pack_index = ', pack_index, 'radiative recombination'
 package(pack_index)%typ = type_kpkt
 count_crecombination = count_crecombination + 1
 active = 0
! no event was chosen
ELSE
 write(*,*) 'do_ipackage, pack_index = ', pack_index, ' no event was chosen...'
END IF

DEALLOCATE(linetransitions, lineuptransitions, &
                Lma_int_dorad, Lma_int_uprad, Lma_int_docoll, Lma_int_upcoll, &
                Lma_int_up, Lma_int_do)
IF(nlevslion /= 0) DEALLOCATE(Lma_recrad, Lma_int_recrad, Lma_reccol, Lma_int_reccol)
END DO

END SUBROUTINE do_ipackage
