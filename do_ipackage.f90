! this subroutine computes an i-packet dynamics
SUBROUTINE do_ipackage(pack_index)

USE TYPES
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index
INTEGER                         :: actual_state
INTEGER                         :: last_line, last_level
! is a macro atom active?
INTEGER                         :: active
! number of line transitions
INTEGER                         :: nlns, nluns
INTEGER                         :: element_index, ion_index
INTEGER                         :: I, J, K, line
INTEGER, ALLOCATABLE            :: linetransitions(:), transitions(:), lineuptransitions(:), &
                                   lineradtransitions(:) ! an array of transitions down from last_line
DOUBLE PRECISION, ALLOCATABLE   :: Lintdownrad(:), Lraddeexc(:), Lintuprad(:), &
                                   Lrad(:), Lintdowncoll(:), Lintupcoll(:), Lcolldeexc(:), &
                                   Lintdown(:), Lintup(:)
INTEGER                         :: act_line
DOUBLE PRECISION                :: actVal
DOUBLE PRECISION                :: ran2, rand
! sum function
DOUBLE PRECISION                :: Z, Ztotal, Zintdownrad, Zraddeexc, Zintuprad, Zrad, Zcoll, &
                                   Zintupcoll, Zintdowncoll, Zdown, Zup, Zintdown, Zintup
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1, Z2, Z3
DOUBLE PRECISION                :: summ, stat_weight, exci_energy
! populations
DOUBLE PRECISION                :: act_pop
INTEGER                         :: get_package_model_index, current_mgi
! new frequency
DOUBLE PRECISION                :: new_freq
! Doppler factor
DOUBLE PRECISION                :: D

! define the needed variables
last_line = package(pack_index)%last_line
last_level = linelist(last_line)%upper
element_index = linelist(last_line)%indexe
ion_index = linelist(last_line)%indexi
current_mgi = get_package_model_index(pack_index)
print*, 'do_ipackage: element_index = ', element_index, 'ion_index = ', ion_index

active = 1
! this is an initial state of the macro-atom
actual_state = last_level
! we will run this loop until the macro atom is deactivated
DO WHILE (active == 1)
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
 ALLOCATE(Lraddeexc(nlns))
 ! 3.) collisional deexcitation
 ALLOCATE(Lintdowncoll(nlns))
 ALLOCATE(Lcolldeexc(nlns))
 ! 2.) internal upward jump within the current ion
 ALLOCATE(Lintuprad(nluns))
 ALLOCATE(Lintupcoll(nluns))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! calculation of the given transition probabilities
 CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, act_pop)
 CALL radiative_rates(nlns, linetransitions, nluns, lineuptransitions, act_pop, &
 Lintdownrad, Zintdownrad, Lintuprad, Zintuprad, Lraddeexc, Zraddeexc)
 CALL collisional_rates(1, pack_index, actual_state, nlns, linetransitions, nluns, lineuptransitions, &
 act_pop, Lintdowncoll, Zintdowncoll, Lintupcoll, Zintupcoll, Lcolldeexc, Zcoll)
 ! total rates of internal donwnward jump
 DO I = 1, nluns
   Lintdown(I) = Lintdownrad(I) + Lintdowncoll(I)
   Zintdown = Zintdownrad + Zintdowncoll
 END DO
 ! total rates of internal upward jump
 DO I = 1, nluns
  ! internal jump up
   Lintup(I) = Lintuprad(I) + Lintupcoll(I)
   Zintup = Zintuprad + Zintupcoll
 END DO
! the total sum 
Ztotal = Zintdown+ Zraddeexc + Zintup + Zcoll
! a random number for computation, which process occurs
rand = ran2(idum)
!print*, 'do_ipackage: Ztotal = ', Ztotal
rand = rand * Ztotal
!print*, 'do_ipackage: random number: ', rand
! these variables are only to the whole line won't be too long
Z0 = Zintdown
Z2 = Zintup
Z1 = Zraddeexc
Z3 = Zcoll
!print*, 'Zintdownjump, Zintupjump, Zraddeexc, Zcoll: ', Zintdownjump, Zintupjump, Zraddeexc, Zcoll
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
   print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
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
 print*, 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
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
  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == ion_index) THEN
   ! transitions to a lower level
   IF(linelist(K)%upper == linelist(last_line)%upper) THEN
    nlns = nlns + 1
   END IF
  END IF
 END DO
 ALLOCATE(lineradtransitions(nlns), Lrad(nlns))
 J = 0
 DO K = 1, ntransitions
  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == ion_index) THEN
   ! transitions to a lower level
   IF(linelist(K)%upper == linelist(last_line)%upper) THEN
    J = J + 1
    lineradtransitions(J) = K
    actVal = linelist(K)%A_ul
    Lrad(J) = actVal
    Zrad = Zrad + actVal
   END IF
  END IF
 END DO
 rand = rand * Zrad
 !print*, 'do_ipackage: rand = ', rand, ' Zrad = ', Zrad
 summ = 0
 ! looking for the given line
 DO line = 1, nlns
  IF(rand >= summ .AND. rand <= summ + Lrad(line)) THEN
!   print*, 'radiative deexcitation'
   ! we found the given cell now we have to compute only a new frequency
   new_freq = linelist(lineradtransitions(line))%freq
   package(pack_index)%freq_cmf = new_freq
   CALL doppler_factor(pack_index, D)
   package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
   package(pack_index)%last_line = no_line
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
   print*, 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
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
 print*, 'pack_index = ', pack_index, ' collisional deexcitation...'
! DO I = 1, nlns
!  ! we will find the given state
!  IF(rand >= summ .AND. rand < summ + Ldowncoll(I)) THEN
!   !print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
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
END IF


! only one loop
!STOP 'testing the code'

DEALLOCATE(linetransitions, lineuptransitions, &
                Lintdownrad, Lintuprad, Lraddeexc, Lintdowncoll, Lintupcoll, &
                Lcolldeexc, Lintup, Lintdown)

END DO

END SUBROUTINE do_ipackage
