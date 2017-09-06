! this subroutine computes an i-packet dynamics
SUBROUTINE do_ipackage(pack_index)

USE TYPES
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
! RATES FOR TRANSITIONS
! Lintdownrad -- internal downward jumps: radiative part
! Lintuprad -- internal upward jumps: radiative part
! Lrad -- radiative transitions, computed if and only if radiative deexcitation occurs
!         downward transitions from the last_line
! Lintdowncoll -- internal downward jumps: collisional part
! Lintupcoll -- internal upward jumps: collisional part
! Lintdown -- internal downward jumps: radiative + collisional
! Lintup -- internal upward jumps: radiative + collisional
! Lphotrecom -- photorecombination
! Lintphotrecom -- internal photorecombination
! Lcollrecom -- collisional recombination
! Lintcollrecom -- internal collisional recombination
DOUBLE PRECISION, ALLOCATABLE   :: Lintdownrad(:), Lintuprad(:), Lrad(:), &
                                   Lintdowncoll(:), Lintupcoll(:), &
                                   Lintdown(:), Lintup(:), &
                                   ! internal recombination
                                   Lphotrecom(:), Lintphotrecom(:), Lcollrecom(:), Lintcollrecom(:)
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
DOUBLE PRECISION                :: Z0, Z1, Z2, Z3, Z4, Z5, Z6
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
 ALLOCATE(Lintdownrad(nlns))
 ! 1.) radiative deexcitation
 ! this allocates only in the first loop, because this field will only remember transitions
 ! from the last_line's upper level
 ! 3.) collisional deexcitation
 ALLOCATE(Lintdowncoll(nlns))
 ! 2.) internal upward jump within the current ion
 ALLOCATE(Lintuprad(nluns))
 ALLOCATE(Lintupcoll(nluns))
 ALLOCATE(Lintdown(nlns), Lintup(nluns))
 Zintdown = 0.D0
 Zintup = 0.D0
 ! allocation of the field for recombination processes
 IF(ion_index > 1) THEN
  nlevslion = SIZE(elements(element_index)%ions(ion_index - 1)%levels)
  !print*, 'element_index = ', element_index, 'ion_index - 1 = ', ion_index - 1, ' nlevslion = ', nlevslion
  ALLOCATE(Lphotrecom(nlevslion), Lintphotrecom(nlevslion), Lcollrecom(nlevslion), Lintcollrecom(nlevslion))
 ELSE
  nlevslion = 0
 END IF
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! calculation of the given transition probabilities
 CALL populations(element_index, ion_index, actual_state, current_mgi, act_pop)
 CALL radiative_rates(nlns, linetransitions, nluns, lineuptransitions, act_pop, &
  Lintdownrad, Zintdownrad, Lintuprad, Zintuprad, Zraddeexc)
 CALL collisional_rates(1, pack_index, actual_state, nlns, linetransitions, nluns, lineuptransitions, &
  act_pop, Lintdowncoll, Zintdowncoll, Lintupcoll, Zintupcoll, Zcoll)
 CALL photion_rates(0, element_index, ion_index, actual_state, current_mgi, act_pop, Zphotionup, &
  nlevslion, Lphotrecom, Zphotrecom, Lintphotrecom, Zphotiondown)
 CALL collion_rates(1, element_index, ion_index, pack_index, actual_state, act_pop, Zcollionup, &
  nlevslion, Lintcollrecom, Zcolliondown,Lcollrecom, Zcollrecom)
 ! total rates of internal donwnward jump
 DO I = 1, nlns
   Lintdown(I) = Lintdownrad(I) + Lintdowncoll(I)
 END DO
 Zintdown = Zintdownrad + Zintdowncoll
 ! total rates of internal upward jump
 DO I = 1, nluns
  ! internal jump up
   Lintup(I) = Lintuprad(I) + Lintupcoll(I)
   !print*, 'Lintuprad(I) = ', Lintuprad(I), ' Lintupcoll(I) = ', Lintupcoll(I)
 END DO
 Zintup = Zintuprad + Zintupcoll
! an internal ionization sum
Zionization = Zphotionup + Zcollionup
Zintrecombination = Zphotiondown + Zcolliondown
Zrecombination = Zphotrecom + Zcollrecom
! the total sum 
Ztotal = Zintdown + Zraddeexc + Zintup + Zcoll + Zionization + Zrecombination + Zintrecombination
! a random number for computation, which process occurs
rand = ran2(idum)
!print*, 'do_ipackage: Ztotal = ', Ztotal
rand = rand * Ztotal
! print*, 'do_ipackage: random number: ', rand
! these variables are only to the whole line won't be too long
Z0 = Zintdown
Z1 = Z0 + Zraddeexc
Z2 = Z1 + Zintup
Z3 = Z2 + Zcoll
Z4 = Z3 + Zionization
Z5 = Z4 + Zintrecombination
Z6 = Z5 + Zrecombination
!print*, 'Zintdown, Zintup, Zraddeexc, Zcoll, Zionization, Zrecombination, Zintrecombination: ', &
!        Zintdown, Zintup, Zraddeexc, Zcoll, Zionization, Zrecombination, Zintrecombination
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump
! in this case a macro-atom transits into a lower state without an energy emission
IF(rand >= 0 .AND. rand < Z0) THEN
 count_intdownjump = count_intdownjump + 1
! next transition will be an internal downward jump
 summ = 0
 DO I = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintdown(I)) THEN
   actual_state = linelist(linetransitions(I))%lower
   !print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
   EXIT
  END IF
  summ = summ + Lintdown(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radiative deexciation
! in this case the packet is now transformed back to the r-packet
! we have to choose a new frequency, which would be calculated randomly from the 
! possible transition last_line -> some lower line
ELSE IF (rand >= Z0 .AND. rand <= Z1) THEN
! next transition will be an radiative deexcitation
 !print*, 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
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
 ALLOCATE(lineradtransitions(nlns), Lrad(nlns))
 J = 0
! print*, 'do_ipackage: nlns = ', nlns, &
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
    Lrad(J) = linelist(K)%A_ul * stat_weight * (exci_energy_u - exci_energy_l)
    !print*, 'do_ipackage: Lrad(J) = ', Lrad(J)
    Zrad = Zrad + Lrad(J)
   END IF
  END IF
 END DO
 rand = rand * Zrad
 summ = 0
 ! looking for the given line
 DO line = 1, nlns
  IF(rand >= summ .AND. rand <= summ + Lrad(line)) THEN
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
  summ = summ + Lrad(line)
 END DO
 active = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF (rand >= Z1 .AND. rand <= Z2) THEN
 count_intupjump = count_intupjump + 1
 summ = Z1
 DO I = 1, nluns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintup(I)) THEN
   actual_state = linelist(lineuptransitions(I))%upper
   !print*, 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
   EXIT
  END IF
  summ = summ + Lintup(I)
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
  !print*, 'summ = ', summ, ' summ + L(I) = ', summ + Lphotrecom(I) + Lcollrecom(I)
  IF(rand >= summ .AND. rand < summ + Lintphotrecom(I) + Lintcollrecom(I)) THEN
   actual_state = I
!   print*, 'do_ipackage: changing actual state to the state I = ', I
   !print*, 'do_ipackage: packet: ', pack_index, ' internal jump to the lower ionization state...'
   EXIT
  END IF
  summ = summ + Lintphotrecom(I) + Lintcollrecom(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! recombination (radiative and collisional)
ELSE IF(rand >= Z5 .AND. rand <= Z6) THEN
 !print*, 'do_ipackage: pack_index = ', pack_index, 'radiative recombination'
 package(pack_index)%typ = type_kpkt
 active = 0
 ! temporary
 ! the frequency should be sampled from the photion cross section
 count_recombination = count_recombination + 1
 CALL i_freq_recomb(actual_state, pack_index, act_pop, new_freq)
 package(pack_index)%freq_cmf = new_freq
END IF

DEALLOCATE(linetransitions, lineuptransitions, &
                Lintdownrad, Lintuprad, Lintdowncoll, Lintupcoll, &
                Lintup, Lintdown)
IF(nlevslion /= 0) DEALLOCATE(Lphotrecom, Lintphotrecom, Lcollrecom, Lintcollrecom)
END DO

END SUBROUTINE do_ipackage
