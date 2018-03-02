! this subroutine computes an i-packet dynamics
SUBROUTINE do_ipackage(pack_index)

USE TYPES
USE rates_i
USE counters
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index
INTEGER                         :: actual_state
INTEGER                         :: last_line, last_level, last_ion
! is a macro atom active?
INTEGER                         :: active
! number of line transitions
INTEGER                         :: nlns, nluns, act_line
INTEGER                         :: element_index, ion_index
INTEGER                         :: I, J, K, line
INTEGER, ALLOCATABLE            :: linetransitions(:), transitions(:), lineuptransitions(:), &
                                   lineradtransitions(:) ! an array of transitions down from last_line
INTEGER                         :: nlevslion
DOUBLE PRECISION                :: actVal
DOUBLE PRECISION                :: rand
! sum function
DOUBLE PRECISION                :: Z, Ztotal, Zintdownrad, Zraddeexc, Zintuprad, Zrad, Zcoll, &
                                   Zintupcoll, Zintdowncoll, Zdown, Zup, Zintdown, Zintup, &
                                   Zphotrecom, Zintrecombination, Zcollrecom
! B-F internal processes
DOUBLE PRECISION                :: Zphotionup, Zphotiondown, Zcollionup, Zcolliondown
DOUBLE PRECISION                :: Zionization, Zrecombination
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7
DOUBLE PRECISION                :: summ, stat_weight_u, stat_weight_l, exci_energy, exci_energy_l, exci_energy_u
! populations
DOUBLE PRECISION                :: act_pop, low_pop
INTEGER                         :: get_package_model_index, current_mgi
! new frequency
DOUBLE PRECISION                :: new_freq
! beta calculation
DOUBLE PRECISION                :: taulu, betalu, Blu
! Doppler factor
DOUBLE PRECISION                :: D
TYPE(irates)                    :: actirates
INTEGER                         :: OMP_GET_THREAD_NUM, my_rank
REAL(8)                         :: random
! write down the processes
LOGICAL                         :: procout = .FALSE.

my_rank = OMP_GET_THREAD_NUM()

last_line = package(pack_index)%last_line
! define the needed variables
! it is necessary to remember the initial conditions of a macro-atom
!IF(.NOT. ASSOCIATED(actirates)) ALLOCATE(actirates)
element_index = package(pack_index)%l_ele
last_ion = package(pack_index)%l_ion
last_level = package(pack_index)%l_lev
IF(package(pack_index)%last_line /= no_line) linelist(last_line)%n_exc = linelist(last_line)%n_exc + 1
 
current_mgi = get_package_model_index(pack_index)
! write(*,*) '*********************************************************************'
! write(*,*)  'do_ipackage: element_index = ', element_index, 'ion_index = ', last_ion, ' level_index = ', last_level

active = 1
! this is an initial state of the macro-atom
actual_state = last_level
ion_index = last_ion
! we will run this loop until the macro atom is deactivated
DO WHILE (active == 1)
 ! write(*,*) 'do_ipackage: actual_state = ', actual_state, ' ion_index = ', ion_index
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
 !ALLOCATE(actirates%Lma_int_dorad(nlns))
 ! 1.) radiative deexcitation
 ! this allocates only in the first loop, because this field will only remember transitions
 ! from the last_line's upper level
 ! 3.) collisional deexcitation
 !ALLOCATE(actirates%Lma_int_docoll(nlns))
 ! 2.) internal upward jump within the current ion
 !ALLOCATE(actirates%Lma_int_uprad(nluns))
 !ALLOCATE(actirates%Lma_int_upcoll(nluns))
 !ALLOCATE(actirates%Lma_int_do(nlns), actirates%Lma_int_up(nluns))
 Zintdown = 0.D0
 Zintup = 0.D0
 !IF(my_rank == 1) THEN
 ! write(*,*) 'do_ipackage: nlns = ', nlns, ' nluns = ', nluns
 ! write(*,*) 'do_ipackage: int_do = ', size(actirates%Lma_int_dorad), &
 !            ' int_docoll = ', size(actirates%Lma_int_docoll), &
 !            ' intuprad = ', size(actirates%Lma_int_uprad), &
 !            ' intupcoll = ', size(actirates%Lma_int_upcoll), &
 !            ' int_do = ', size(actirates%Lma_int_do), &
 !            ' intup = ', size(actirates%Lma_int_up)
 !END IF
 ! allocation of the field for recombination processes
 IF(ion_index > 1) THEN
  nlevslion = SIZE(elements(element_index)%ions(ion_index - 1)%levels)
 ! write(*,*) 'do_ipackage: indexe = ', element_index, 'indexi - 1 = ', &
 !  ion_index - 1, ' nlevslion = ', nlevslion
 ELSE
  nlevslion = 0
 END IF
 actirates = irates(nlns, nluns, nlevslion)
 !write(*,*) 'ALLOCATED: do_ipackage: my_rank = ', my_rank, ' recrad = ', size(actirates%Lma_recrad), &
 !           ' intdorad = ', size(actirates%Lma_int_dorad)
 !write(*,*) 'do_ipackage: loc(actirates) = ', loc(actirates)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! calculation of the given transition probabilities
  ! write(*,*) 'do_ipackage: element_index = ', element_index, ' ion_index = ', ion_index
 CALL populations(element_index, ion_index, actual_state, current_mgi, act_pop)
 ! write(*,*) 'do_ipackage: pop = ', act_pop
 CALL i_radtrans(current_mgi, nlns, linetransitions, nluns, lineuptransitions, act_pop, &
  Zintdownrad, Zintuprad, Zraddeexc, actirates)
 CALL i_coltrans(1, pack_index, actual_state, nlns, linetransitions, nluns, lineuptransitions, &
  act_pop, Zintdowncoll, Zintupcoll, Zcoll, actirates)
 CALL i_radion(0, element_index, ion_index, actual_state, current_mgi, act_pop, Zphotionup, &
   Zphotiondown, Zphotrecom, actirates)
 CALL i_colion(1, element_index, ion_index, actual_state, pack_index, act_pop, Zcollionup, &
  Zcolliondown,Zcollrecom, actirates)
! testing
! Zraddeexc = 0.D0
! controll part
IF(Zintdowncoll < 0.D0) STOP 'do_ipackage: Zintdowncoll < 0'
IF(Zintupcoll < 0.D0) STOP 'do_ipackage: Zintupcoll < 0'
IF(Zcoll < 0.D0) STOP 'do_ipackage:  Zcoll < 0'
IF(Zphotionup < 0.D0) STOP 'do_ipackage: Zphotionup < 0'
IF(Zphotiondown < 0.D0) STOP 'do_ipackage: Zphotiondown < 0'
IF(Zphotrecom < 0.D0) STOP 'do_ipackage: Zphotrecom < 0'
IF(Zcolliondown < 0.D0) STOP 'do_ipackage: Zcolliondown < 0'
IF(Zcollrecom < 0.D0) STOP 'do_ipackage: Zcollrecom < 0'
 ! total rates of internal donwnward jump
 DO I = 1, nlns
   actirates%Lma_int_do(I) = actirates%Lma_int_dorad(I) + actirates%Lma_int_docoll(I)
 END DO
 Zintdown = Zintdownrad + Zintdowncoll
 ! total rates of internal upward jump
 DO I = 1, nluns
  ! internal jump up
   actirates%Lma_int_up(I) = actirates%Lma_int_uprad(I) + actirates%Lma_int_upcoll(I)
   ! write(*,*)  'actirates%Lma_int_uprad(I) = ', actirates%Lma_int_uprad(I), &
   !  ' actirates%Lma_int_upcoll(I) = ', actirates%Lma_int_upcoll(I)
 END DO
 Zintup = Zintuprad + Zintupcoll
! an internal ionization sum
Zionization = Zphotionup + Zcollionup
 ! write(*,*) 'do_ipackage: Zphotionup = ', Zphotionup, ' Zcollionup = ', Zcollionup
Zintrecombination = Zphotiondown + Zcolliondown
Zrecombination = Zphotrecom + Zcollrecom
! write(*,*) 'do_ipackage: Zphotrecom = ', Zphotrecom
! the total sum 
Ztotal = Zintdown + Zraddeexc + Zintup + Zcoll + &
    Zionization + Zrecombination + Zintrecombination
! a random number for computation, which process occurs
rand = DBLE(random()) * Ztotal
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
! write(*,*) 'Zintdown, Zintup, Zraddeexc, Zcoll, Zionization, Zrecombination, Zintrecombination, Zphotrecom: ', &
!         Zintdown, Zintup, Zraddeexc, Zcoll, Zionization, Zrecombination, Zintrecombination, Zphotrecom
!write(*,*) 'do_ipackage: Zrecombination = ', Zrecombination, ' Zphotrecom = ', Zphotrecom, ' Zcollrecom = ', Zcollrecom
! write(*,*) 'do_ipackage: Z0 = ', Z0, ' Z1 = ', Z1, ' Z2 = ', Z2, ' Z3 = ', Z3, &
!  ' Z4 = ', Z4, ' Z5 = ', Z5, ' Z6 =', Z6, ' Z7 = ', Z7
IF(Ztotal == 0.D0) STOP 'do_ipackage: Ztotal = 0'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump
! in this case a macro-atom transits into a lower state without an energy emission
IF(rand >= 0.D0 .AND. rand < Z0) THEN
 !$OMP ATOMIC
 count_i_int_down = count_i_int_down + 1
! next transition will be an internal downward jump
 summ = 0.D0
 DO I = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + actirates%Lma_int_do(I)) THEN
   actual_state = linelist(linetransitions(I))%lower
   IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
   EXIT
  END IF
  summ = summ + actirates%Lma_int_do(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radiative deexciation
! in this case the packet is now transformed back to the r-packet
! we have to choose a new frequency, which would be calculated randomly from the 
! possible transition last_line -> some lower line
ELSE IF (rand >= Z0 .AND. rand <= Z1) THEN
! next transition will be an radiative deexcitation
 package(pack_index)%typ = type_rpkt
 ! IF(package(pack_index)%last_line == no_line) CYCLE
 package(pack_index)%last_line = no_line
 CALL emit_rpackage(pack_index)
 IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
 ! now we will calculate new frequency of the packet
 ! we will choose this frequency from the possible radiative transitions
 rand = random()
 Zrad = 0.D0
 nlns = 0
! DO K = 1, ntransitions
!  ! we are interested only in the transitions for the given atom
!  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == ion_index) THEN
!   ! transitions to a lower level
!  ! print*, 'do_ipackage: linelist(K)%upper = ', linelist(K)%upper
!   IF(linelist(K)%upper == last_level) THEN
!    nlns = nlns + 1
!   END IF
!  END IF
! END DO
! ALLOCATE(lineradtransitions(nlns), actirates%Lma_rad(nlns))
! ALLOCATE(actirates%Lma_rad(nlns))
! J = 0
!  ' linelist(last_line)%indexe = ', linelist(last_line)%indexe,&
!  ' linelist(last_line)%indexi = ', linelist(last_line)%indexi, ' last_level = ', last_level, &
!  ' last_line = ', last_line, ' ntransitions = ', ntransitions
! DO K = 1, ntransitions
!  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == ion_index) THEN
!   ! transitions to a lower level
!   IF(linelist(K)%upper == last_level) THEN
!    J = J + 1
! DO K = 1, nlns
!  ! lineradtransitions(J) = K
!  act_line = linetransitions(K)
!  exci_energy_u = &
!   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
!  exci_energy_l = &
!   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
!  stat_weight_l = &
!   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%stat_waight
!  stat_weight_u = &
!   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
!  CALL populations(element_index, ion_index, linelist(act_line)%lower, current_mgi, low_pop)
!  Blu = 4.0 * pi / (h * linelist(K)%freq) * linelist(act_line)%A_ul
!  taulu = low_pop * Blu * h * light_speed / (4.0 * pi) *&
!   (1.D0 - (stat_weight_l * act_pop) / (stat_weight_u * low_pop))
!  betalu = 1 / taulu * (1 - exp(-taulu))
!  actirates%Lma_rad(K)  = act_pop * betalu * linelist(K)%A_ul * &
!   (exci_energy_u - exci_energy_l)
!  !print*, 'do_ipackage: Lrad(J) = ', Lrad(J)
!  Zrad = Zrad + actirates%Lma_rad(K)
!  ! END IF
!  !END IF
! END DO
! rand = rand * Zrad
! summ = 0
 ! looking for the given line
 DO line = 1, nlns
  IF(rand >= summ .AND. rand <= summ + actirates%Lma_rad(line)) THEN
   IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
   !print*, 'raddeexc: summ = ', summ, ' rand = ', rand, ' Zrad = ', Zrad
   ! we found the given cell now we have to compute only a new frequency
   new_freq = linelist(linetransitions(line))%freq
   CALL emit_rpackage(pack_index)
   package(pack_index)%freq_cmf = new_freq
   CALL doppler_factor(pack_index, D)
   package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
   package(pack_index)%last_line = no_line
   !$OMP ATOMIC
   linelist(linetransitions(line))%n_deexc = linelist(linetransitions(line))%n_deexc + 1
   IF(last_line /= no_line) THEN
    IF(linelist(linetransitions(line))%lower == linelist(last_line)%lower) THEN
     ! resonant scattering occures
     !$OMP ATOMIC
     count_i_rad_dxrs = count_i_rad_dxrs + 1
    ELSE
     !$OMP ATOMIC
     count_i_rad_dxfl = count_i_rad_dxfl + 1
    END IF
   END IF
  END IF
    EXIT
  summ = summ + actirates%Lma_rad(line)
 END DO
 active = 0
 DEALLOCATE(actirates%Lma_rad)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF (rand >= Z1 .AND. rand <= Z2) THEN
 !$OMP ATOMIC
 count_i_int_upwa = count_i_int_upwa + 1
 summ = Z1
 DO I = 1, nluns
  ! we will find the given state
  !write(*,*) 'do_ipackage: summ = ', summ, ' rand = ', rand, ' summ + actirates%Lma_int_up = ', summ + actirates%Lma_int_up(I)
  IF(rand >= summ .AND. rand < summ + actirates%Lma_int_up(I)) THEN
   actual_state = linelist(lineuptransitions(I))%upper
   IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
   EXIT
  END IF
  summ = summ + actirates%Lma_int_up(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional deexcitation
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF(rand >= Z2 .AND. rand <= Z3) THEN
 summ = Z2
 package(pack_index)%last_line = no_line
 package(pack_index)%typ = type_kpkt
 IF(procout) write(*,*)  'pack_index = ', pack_index, ' collisional deexcitation...'
 active = 0
 !$OMP ATOMIC
 count_i_col_deex = count_i_col_deex + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal photoionization
ELSE IF(rand >= Z3 .AND. rand <= Z4) THEN
 IF(procout) write(*,*)  'pack_index = ', pack_index, ' internal jump to to the upper ionization state...'
 ion_index = ion_index + 1
 !$OMP ATOMIC
 count_i_int_phot = count_i_int_phot + 1
 actual_state = 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal recombination
ELSE IF(rand >= Z4 .AND. rand <= Z5) THEN
 IF(procout) write(*,*)  'pack_index = ', pack_index, ' internal jump to to the lower ionization state...'
 ion_index = ion_index - 1
 summ = Z4
 DO I = 1, nlevslion
  ! we will find the given state
  ! print*, 'summ = ', summ, ' summ + L(I) = ', summ + actirates%Lma_int_recrad(I) + actirates%Lma_int_reccol(I)
  IF(rand >= summ .AND. rand < summ + actirates%Lma_int_recrad(I) + actirates%Lma_int_reccol(I)) THEN
   actual_state = I
   !$OMP ATOMIC
   count_i_int_reco = count_i_int_reco + 1
   IF(procout) write(*,*)  'do_ipackage: packet: ', pack_index, ' internal jump to the lower ionization state...', 'actual_state = &
   ', actual_state
   EXIT
  END IF
  summ = summ + actirates%Lma_int_recrad(I) + actirates%Lma_int_reccol(I)
 END DO
 ! write(*,*) 'do_ipackage: Z5 = ', Z5, ' summ = ', summ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radiative recombination
ELSE IF(rand >= Z5 .AND. rand <= Z6) THEN
 package(pack_index)%typ = type_rpkt
 active = 0
 ! the frequency should be sampled from the photion cross section
 summ = Z5
 DO I = 1, nlevslion
  ! write(*,*) 'do_ipackage: summ = ', summ, ' rand = ', rand, ' summ + Lma_recrad = ', summ + actirates%Lma_recrad(I)
  IF( rand >= summ .AND. rand < summ + actirates%Lma_recrad(I)) THEN
   IF(procout) write(*,*) 'do_ipackage: pack_index = ', pack_index, 'radiative recombination'
   CALL i_freq_recomb(element_index, ion_index, I, pack_index, act_pop, new_freq)
   package(pack_index)%freq_cmf = new_freq
   CALL doppler_factor(pack_index, D)
   CALL emit_rpackage(pack_index)
   package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
   !$OMP ATOMIC
   count_i_rad_reco = count_i_rad_reco + 1
   EXIT
  END IF
  summ = summ + actirates%Lma_recrad(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional recombination
ELSE IF(rand >= Z6 .AND. rand <= Z7) THEN
 IF(procout) write(*,*) 'do_ipackage: pack_index = ', pack_index, 'collisional recombination'
 package(pack_index)%typ = type_kpkt
 active = 0
 !$OMP ATOMIC
 count_i_col_reco = count_i_col_reco + 1
! no event was chosen
ELSE
 write(*,*) 'do_ipackage, pack_index = ', pack_index, ' no event was chosen...'
 STOP
END IF

DEALLOCATE(linetransitions, lineuptransitions, &
                actirates%Lma_int_dorad, actirates%Lma_int_uprad, actirates%Lma_int_docoll, actirates%Lma_int_upcoll, &
                actirates%Lma_int_up, actirates%Lma_int_do)
IF(nlevslion /= 0) DEALLOCATE(actirates%Lma_recrad, actirates%Lma_int_recrad, actirates%Lma_reccol, actirates%Lma_int_reccol)
END DO

! deallocate rates
IF(ALLOCATED(actirates%Lma_int_dorad)) DEALLOCATE(actirates%Lma_int_dorad)
IF(ALLOCATED(actirates%Lma_int_uprad)) DEALLOCATE(actirates%Lma_int_uprad)
IF(ALLOCATED(actirates%Lma_int_dorad)) DEALLOCATE(actirates%Lma_rad)
IF(ALLOCATED(actirates%Lma_int_docoll)) DEALLOCATE(actirates%Lma_int_docoll)
IF(ALLOCATED(actirates%Lma_int_upcoll)) DEALLOCATE(actirates%Lma_int_upcoll)
IF(ALLOCATED(actirates%Lma_int_up)) DEALLOCATE(actirates%Lma_int_up)
IF(ALLOCATED(actirates%Lma_int_do)) DEALLOCATE(actirates%Lma_int_do)
IF(ALLOCATED(actirates%Lma_recrad)) DEALLOCATE(actirates%Lma_recrad)
!write(*,*) 'do_ipackage: size1 = ', SIZE(actirates%Lma_int_recrad)
IF(ALLOCATED(actirates%Lma_int_recrad)) DEALLOCATE(actirates%Lma_int_recrad)
!write(*,*) 'do_ipackage: size2 = ', SIZE(actirates%Lma_int_reccol)
IF(ALLOCATED(actirates%Lma_int_reccol)) DEALLOCATE(actirates%Lma_int_reccol)
IF(ALLOCATED(actirates%Lma_int_dorad)) DEALLOCATE(actirates%Lma_int_reccol)


END SUBROUTINE do_ipackage
