! this subroutine computes an i-packet dynamics
! it follows the macroatom approach Lucy(2002, 2003)
! possible transitions:
!  1. internal downvard jump i-pack -> i-pack
!  2. radiative deexcitation i-pack -> r-pack
!  3. internal upvard jump i-pack -> i-pack
!  4. collisional deexcitation i-pack -> k-pack
!  5. internal ionization i-pack -> i-pack
!  6. internal recombination i-pack -> i-pack
!  7. radiative recombination i-pack -> r-pack
!  8. collisional recombination i-pack -> k-pack
!
! INPUT: pack_index(INT): package index
! OUTPUT: NONE
!
SUBROUTINE do_ipackage(pack_index)

USE TYPES
USE constants
USE rates_i
USE counters
USE dummypacket
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
INTEGER                         :: ind_I, line
INTEGER, ALLOCATABLE            :: linetransitions(:), lineuptransitions(:)
INTEGER                         :: nlevslion
DOUBLE PRECISION                :: rand
! sum function
DOUBLE PRECISION                :: Ztotal, Zintdownrad, Zraddeexc, Zintuprad, Zcoll, &
                                   Zintupcoll, Zintdowncoll, Zintdown, Zintup, &
                                   Zphotrecom, Zintrecombination, Zcollrecom
! B-F internal processes
DOUBLE PRECISION                :: Zphotionup, Zphotiondown, Zcollionup, Zcolliondown
DOUBLE PRECISION                :: Zionization, Zrecombination
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1, Z2, Z3, Z4, Z5, Z6, Z7
DOUBLE PRECISION                :: summ
DOUBLE PRECISION                :: ran2
! populations
DOUBLE PRECISION                :: act_pop
INTEGER                         :: get_package_model_index, current_mgi
! new frequency
DOUBLE PRECISION                :: new_freq
! Doppler factor
DOUBLE PRECISION                :: doppler_D
TYPE(irates)                    :: actirates
! write down the processes
LOGICAL                         :: procout = .FALSE.
LOGICAL                         :: sstates = .FALSE.
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: new_dir
INTEGER                         :: dummy_pack_index
! number of processes in MA
! INTEGER, PARAMETER              :: maxproc = 1000000
! INTEGER                         :: n_proc

IF(debug == 4) procout = .TRUE.

IF(pack_index <= SIZE(package)) THEN
 last_line = package(pack_index)%last_line
 element_index = package(pack_index)%l_ele
 last_ion = package(pack_index)%l_ion
 last_level = package(pack_index)%l_lev
 package(pack_index)%l_ele = 0
 package(pack_index)%l_ion = 0
 package(pack_index)%l_lev = 0
 IF(package(pack_index)%last_line /= no_line) &
  linelist(last_line)%n_exc = linelist(last_line)%n_exc + 1
 package(pack_index)%n_interactions = package(pack_index)%n_interactions + 1
ELSE IF(pack_index > SIZE(package)) THEN
 dummy_pack_index = pack_index - SIZE(package)
 last_line = dummypackage(dummy_pack_index)%last_line
 element_index = dummypackage(pack_index)%l_ele
 last_ion = dummypackage(pack_index)%l_ion
 last_level = dummypackage(pack_index)%l_lev
 dummypackage(dummy_pack_index)%l_ele = 0
 dummypackage(dummy_pack_index)%l_ion = 0
 dummypackage(dummy_pack_index)%l_lev = 0
 dummypackage(dummy_pack_index)%n_interactions = &
  dummypackage(dummy_pack_index)%n_interactions + 1
END IF

! define the needed variables
! it is necessary to remember the initial conditions of a macro-atom
!IF(.NOT. ASSOCIATED(actirates)) ALLOCATE(actirates)
 
current_mgi = get_package_model_index(pack_index)

active = 1
! this is an initial state of the macro-atom
actual_state = last_level
ion_index = last_ion
! we will run this loop until the macro atom is deactivated
DO WHILE (active == 1)

 !____________________________________________________________________________________
 ! we have to find all possible downward upward transitions
 ! firstly we calculate number of these possible transitions
 ! number of transitions to a lower level
 IF(element_index == 0) THEN
  write(*,*) 'do_ipackage: element_index == 0, pack_index = ', pack_index
 END IF
 nlns = &
  SIZE(elements(element_index)%ions(ion_index)%levels(actual_state)%linetransitions)
 nluns = &
  SIZE(elements(element_index)%ions(ion_index)%levels(actual_state)%lineuptransitions)
 ! write(*,*) 'do_ipackage: nlns = ', nlns, ' nluns = ', nluns
 ALLOCATE(linetransitions(nlns), lineuptransitions(nluns))
 linetransitions = &
  elements(element_index)%ions(ion_index)%levels(actual_state)%linetransitions
 lineuptransitions = &
  elements(element_index)%ions(ion_index)%levels(actual_state)%lineuptransitions

 !____________________________________________________________________________________
 ! total rates of procedure
 Zintdown = 0.D0
 Zintup = 0.D0

 ! allocation of the field for recombination processes
 IF(ion_index > 1) THEN
  nlevslion = SIZE(elements(element_index)%ions(ion_index - 1)%levels)
 ELSE
  nlevslion = 0
 END IF
 actirates = irates(nlns, nluns, nlevslion)



 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! calculation of the given transition probabilities
  ! write(*,*) 'do_ipackage: element_index = ', element_index, ' ion_index = ', ion_index
  IF(sstates) write(36,*) 'element_index = ', element_index, ' ion_index = ', ion_index, ' actual_state = ',&
   actual_state
 ! write(*,*) 'do_ipackage: pop = ', act_pop
 ! calculations of radiative rates
 CALL i_radtrans(current_mgi, element_index, ion_index, actual_state, &
  Zintdownrad, Zintuprad, Zraddeexc, actirates, pack_index, new_dir)
 CALL i_coltrans(1, pack_index, element_index, ion_index, actual_state, act_pop, &
  Zintdowncoll, Zintupcoll, Zcoll, actirates)
 CALL i_radion(element_index, ion_index, actual_state, current_mgi, act_pop, Zphotionup, &
   Zphotiondown, Zphotrecom, actirates)
 CALL i_colion(1, element_index, ion_index, actual_state, pack_index, act_pop, Zcollionup, &
  Zcolliondown,Zcollrecom, actirates)
! 
! testing
Zcollrecom = 0.D0
Zcoll = 0.D0
Zphotrecom = 0.D0
Zintupcoll = 0.D0
Zintdowncoll = 0.D0
! Zintdownrad = 0.D0
! Zintuprad = 0.D0
Zphotionup = 0.D0
Zcollionup = 0.D0
! end testing
!
! control part of calculation
IF(Zintdowncoll < 0.D0) STOP 'do_ipackage: Zintdowncoll < 0'
IF(Zintupcoll < 0.D0) STOP 'do_ipackage: Zintupcoll < 0'
IF(Zcoll < 0.D0) STOP 'do_ipackage:  Zcoll < 0'
IF(Zphotionup < 0.D0) STOP 'do_ipackage: Zphotionup < 0'
IF(Zphotiondown < 0.D0) STOP 'do_ipackage: Zphotiondown < 0'
IF(Zphotrecom < 0.D0) STOP 'do_ipackage: Zphotrecom < 0'
IF(Zcolliondown < 0.D0) STOP 'do_ipackage: Zcolliondown < 0'
IF(Zcollrecom < 0.D0) STOP 'do_ipackage: Zcollrecom < 0'


 ! total rates of internal donwnward jump
 DO ind_I = 1, nlns
   actirates%Lma_int_do(ind_I) = actirates%Lma_int_dorad(ind_I) + actirates%Lma_int_docoll(ind_I)
 END DO
 Zintdown = Zintdownrad + Zintdowncoll


 ! total rates of internal upward jump
 DO ind_I = 1, nluns
  ! internal jump up
   actirates%Lma_int_up(ind_I) = actirates%Lma_int_uprad(ind_I) + actirates%Lma_int_upcoll(ind_I)
   ! write(*,*)  'actirates%Lma_int_uprad(ind_I) = ', actirates%Lma_int_uprad(ind_I), &
   !  ' actirates%Lma_int_upcoll(ind_I) = ', actirates%Lma_int_upcoll(ind_I)
 END DO

 Zintup = Zintuprad + Zintupcoll
 ! an internal ionization sum
 Zionization = Zphotionup + Zcollionup
 Zintrecombination = Zphotiondown + Zcolliondown
 Zrecombination = Zphotrecom + Zcollrecom
 ! the total sum 
 Ztotal = Zintdown + Zraddeexc + Zintup + Zcoll + &
     Zionization + Zrecombination + Zintrecombination
 ! a random number for computation, which process occurs
 rand = ran2(idum) * Ztotal




 ! these variables are only to the whole line won't be too long
 Z0 = Zintdown
 Z1 = Z0 + Zraddeexc
 Z2 = Z1 + Zintup
 Z3 = Z2 + Zcoll
 Z4 = Z3 + Zionization
 Z5 = Z4 + Zintrecombination
 Z6 = Z5 + Zphotrecom
 Z7 = Z6 + Zcollrecom

! write(*,*) 'Zintdown = ', Zintdown, ' Zraddeexc = ', Zraddeexc, ' Zintup = ', Zintup, ' Zcoll = ', Zcoll, &
!  ' Zionization = ', Zionization, ' Zintrecombination = ', Zintrecombination, ' Zphotrecom = ', Zphotrecom, &
!  ' Zcollrecom = ', Zcollrecom



 IF(Ztotal == 0.D0) THEN
  write(*,*) 'Zintuprad = ', Zintuprad, ' Zintdownrad = ', Zintdownrad
  write(*,*) 'Lma_int_do = ', actirates%Lma_int_do
  write(*,*) 'Lma_int_up = ', actirates%Lma_int_up
  write(*,*) 'actual_state = ', actual_state, ' ion_index = ', ion_index, ' element_index = ', element_index
  write(*,*) 'do_ipackage: Ztotal = 0'
  CALL abort()
 END IF




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! internal downward jump
 ! in this case a macro-atom transits into a lower state without an energy emission
 IF(rand >= 0.D0 .AND. rand < Z0) THEN
  count_i_int_down = count_i_int_down + 1
  IF(procout) write(*,*) 'do_ipackage: internal downward jump...'
 ! next transition will be an internal downward jump
  summ = 0.D0
  DO ind_I = 1, nlns
   ! write(*,*) 'do_ipackage: summ = ', summ, ' rand = ', rand, ' summ + actirates%Lma_int_do = ', summ + actirates%Lma_int_do(I)
   ! we will find the given state
   IF(rand >= summ .AND. rand < summ + actirates%Lma_int_do(ind_I)) THEN
    actual_state = linelist(linetransitions(ind_I))%lower
    IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
    IF(sstates) write(36, *) 'IDJ ->', actual_state
    EXIT
   END IF
   summ = summ + actirates%Lma_int_do(ind_I)
  END DO
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! radiative deexciation
 ! in this case the packet is now transformed back to the r-packet
 ! we have to choose a new frequency, which would be calculated randomly from the 
 ! possible transition last_line -> some lower line
 ELSE IF (rand >= Z0 .AND. rand <= Z1) THEN
 ! 1. radiative deexcitation
  package(pack_index)%typ = type_rpkt
  package(pack_index)%dir = new_dir
  package(pack_index)%next_cross = NONE
  ! IF(package(pack_index)%last_line == no_line) CYCLE
  ! now we will calculate new frequency of the packet
  ! we will choose this frequency from the possible radiative transitions
  summ = 0.D0
  rand = ran2(idum) * Zraddeexc
  ! write(*,*) 'do_ipackage: lines = ', linetransitions(:)
  ! write(*,*) 'do_ipackage: radRate = ', actirates%Lma_rad(:)
  IF(procout) write(*,*) 'do_ipackage: radiative deexcitation...'
  ! looking for the given line
  DO line = 1, nlns
   ! print*, 'raddeexc: summ = ', summ, ' rand = ', rand, ' Zrad = ', Zraddeexc
   IF(rand >= summ .AND. rand <= summ + actirates%Lma_rad(line)) THEN
   ! write(*,*) 'do_ipackage: summ = ', summ, ' rand = ', rand, ' summ + act = ', summ + actirates%Lma_rad(line)
   IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
    ! write(*,*) 'do_ipackage: wale = ', 1.D8 * const_c / linelist(linetransitions(line))%freq
    ! write(*,*) 'do_ipackage: I = ', I
    ! we found the given cell now we have to compute only a new frequency
    new_freq = linelist(linetransitions(line))%freq
    ! testing
    ! new_freq = const_c / (4.D3 * 1.D-8)
    package(pack_index)%freq_cmf = new_freq
    CALL doppler_factor(pack_index, doppler_D)
    IF(sstates) write(36, *) 'RDEEX linewl = ', 1.D8 * const_c / linelist(linetransitions(line))%freq
    package(pack_index)%freq_rf = package(pack_index)%freq_cmf / doppler_D
    package(pack_index)%e_rf = package(pack_index)%e_cmf / doppler_D
    ! write(37,*) 1.D8 * const_c / package(pack_index)%freq_rf
    ! save the emitted frequency
    linelist(linetransitions(line))%n_deexc = linelist(linetransitions(line))%n_deexc + 1
    count_i_rad_deex = count_i_rad_deex + 1
    IF(last_line /= no_line) THEN
     IF(linetransitions(line) == last_line) THEN
      ! resonant scattering occures
      count_i_rad_dxrs = count_i_rad_dxrs + 1
     ELSE
      count_i_rad_dxfl = count_i_rad_dxfl + 1
     END IF
    END IF
    package(pack_index)%last_line = linetransitions(line)
    ! package(pack_index)%last_line = no_line
    EXIT
   END IF
   summ = summ + actirates%Lma_rad(line)
  END DO
  active = 0
  ! STOP 'do_ipackage: testing'
  DEALLOCATE(actirates%Lma_rad)
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! internal upward jump
 ! in this case a macro-atom transits into a upper state without an energy emission
 ELSE IF (rand >= Z1 .AND. rand <= Z2) THEN
  ! write(*,*) 'do_ipackage: Z1 = ', Z1, ' rand = ', rand, ' Z2 = ', Z2
  count_i_int_upwa = count_i_int_upwa + 1
  summ = Z1
  IF(procout) write(*,*) 'do_ipackage: internal upward jump...'
  DO ind_I = 1, nluns
   ! we will find the given state
   ! write(*,*) 'do_ipackage: summ = ', summ, ' rand = ', rand, ' summ + actirates%Lma_int_up = ', summ + actirates%Lma_int_up(I)
   IF(rand >= summ .AND. rand < summ + actirates%Lma_int_up(ind_I)) THEN
    actual_state = linelist(lineuptransitions(ind_I))%upper
    IF(procout) write(*,*) 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
    IF(sstates) write(36, *) 'IUJ ->', actual_state
    EXIT
   END IF
   summ = summ + actirates%Lma_int_up(ind_I)
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
  count_i_col_deex = count_i_col_deex + 1
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! internal photoionization
 ELSE IF(rand >= Z3 .AND. rand <= Z4) THEN
  IF(procout) write(*,*)  'pack_index = ', pack_index, ' internal jump to to the upper ionization state...'
  ion_index = ion_index + 1
  IF(sstates) write(36, *) 'IPHO ->', actual_state
  count_i_int_phot = count_i_int_phot + 1
  actual_state = 1
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! internal recombination
 ELSE IF(rand >= Z4 .AND. rand <= Z5) THEN
  IF(procout) write(*,*)  'pack_index = ', pack_index, ' internal jump to to the lower ionization state...'
  ion_index = ion_index - 1
  summ = Z4
  DO ind_I = 1, nlevslion
   ! we will find the given state
   IF(rand >= summ .AND. rand < summ + actirates%Lma_int_recrad(ind_I) + actirates%Lma_int_reccol(ind_I)) THEN
    actual_state = ind_I
    count_i_int_reco = count_i_int_reco + 1
    IF(sstates) write(36, *) 'IREC ->', actual_state
    IF(procout) write(*,*)  'do_ipackage: packet: ', pack_index, ' internal jump to the lower ionization state...',&
    'actual_state = ', actual_state
    EXIT
   END IF
   summ = summ + actirates%Lma_int_recrad(ind_I) + actirates%Lma_int_reccol(ind_I)
  END DO
  ! write(*,*) 'do_ipackage: Z5 = ', Z5, ' summ = ', summ
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! radiative recombination
 ELSE IF(rand >= Z5 .AND. rand <= Z6) THEN
  package(pack_index)%typ = type_rpkt
  active = 0
  IF(procout) write(*,*) 'do_ipackage: radiative recombination'
  ! the frequency should be sampled from the photion cross section
  summ = Z5
  DO ind_I = 1, nlevslion
   IF( rand >= summ .AND. rand < summ + actirates%Lma_recrad(ind_I)) THEN
    IF(procout) write(*,*) 'do_ipackage: pack_index = ', pack_index, 'radiative recombination'
    package(pack_index)%last_line = no_line
    CALL i_freq_recomb(element_index, ion_index, ind_I, pack_index, new_freq)
    package(pack_index)%freq_cmf = new_freq
    CALL emit_rpackage(pack_index)
    count_i_rad_reco = count_i_rad_reco + 1
    EXIT
   END IF
   summ = summ + actirates%Lma_recrad(ind_I)
  END DO
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! collisional recombination
 ELSE IF(rand >= Z6 .AND. rand <= Z7) THEN
  IF(procout) write(*,*) 'do_ipackage: pack_index = ', pack_index, 'collisional recombination'
  package(pack_index)%typ = type_kpkt
  active = 0
  count_i_col_reco = count_i_col_reco + 1
 ! no event was chosen
 ELSE
 write(*,*) 'do_ipackage: Zintdown = ', Zintdown, ' Zraddeexc = ', Zraddeexc, ' Zintup = ', Zintup, ' Zcoll = ', Zcoll, &
  ' Zionization = ', Zionization, ' Zintrecombination = ', Zintrecombination, ' Zphotrecom = ', Zphotrecom, &
  ' Zcollrecom = ', Zcollrecom
  write(*,*) 'do_ipackage, pack_index = ', pack_index, ' no event was chosen...'
  STOP
 END IF
 
 DEALLOCATE(actirates%Lma_int_dorad, actirates%Lma_int_uprad, actirates%Lma_int_docoll, actirates%Lma_int_upcoll, &
                 actirates%Lma_int_up, actirates%Lma_int_do)
 IF(nlevslion /= 0) DEALLOCATE(actirates%Lma_recrad, actirates%Lma_int_recrad, actirates%Lma_reccol, actirates%Lma_int_reccol)
 DEALLOCATE(linetransitions, lineuptransitions)
 ! n_proc = n_proc + 1
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

! write(*,*) 'do_ipackage: e_rf = ', package(pack_index)%e_rf

END SUBROUTINE do_ipackage
