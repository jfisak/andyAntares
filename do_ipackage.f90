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
                                   lineradtransitions(:)
DOUBLE PRECISION, ALLOCATABLE   :: Lintdownjump(:), Lraddeexc(:), Lintupjump(:), Lrad(:), Lcoll(:)
INTEGER                         :: act_line
DOUBLE PRECISION                :: actVal
DOUBLE PRECISION                :: ran2, rand
! sum function
DOUBLE PRECISION                :: Z, Ztotal, Zintdownjump, Zraddeexc, Zintupjump, Zrad, Zcoll
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1, Z2
DOUBLE PRECISION                :: summ, stat_weight, exci_energy
! populations
DOUBLE PRECISION                :: act_popup, act_popdown
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
 Zintdownjump = 0.D0
 ALLOCATE(Lintdownjump(nlns))
 ! 1.) radiative deexcitation
 Zraddeexc = 0.D0
 ALLOCATE(Lraddeexc(nlns))
 ! 2.) internal upward jump within the current ion
 Zintupjump = 0.D0
 ALLOCATE(Lintupjump(nluns))
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! internal downward jump and radiative deexcitation
 DO I = 1, nlns
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! internal downward jump
  act_line = linetransitions(I)
  ! internal jump down
  ! calculation of number density of the given ion
  CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, act_popdown)
  ! statistical weight
  stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
  exci_energy = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
  ! calculation of a rate coefficient
  actVal = act_popdown * linelist(act_line)%A_ul * exci_energy
 ! print*, 'do_ipackage: stat_waight, exci_energy, act_popdown, linelist(act_line)%A_ul, actVal', &
 !       stat_weight, exci_energy, act_popdown, linelist(act_line)%A_ul, actVal
  Zintdownjump = Zintdownjump + stat_weight * actVal
  Lintdownjump(I) = actVal
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! radiative deexcitation
  CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, act_popup)
  !print*, 'do_ipackage: act_popup = ', act_popup
  actVal = act_popup * linelist(act_line)%A_ul * &
   (elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy - &
   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy)
  Zraddeexc = Zraddeexc + stat_weight * actVal
!  print*, 'do_ipackage: e_u - e_l, actVal', &
!   (elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy - &
!   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy), &
!   actVal
  Lraddeexc(I) = actVal
 END DO
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! collisional deexcitacion
 ! now in the van Regemorter approximation (the first parameter is equal to 1)
 ALLOCATE(Lcoll(nlns))
 CALL collisional_rates(1, pack_index, act_line, nlns, linetransitions, Zcoll, Lcoll)
!print*, 'do_ipackage: Zintdownjump = ', Zintdownjump
 DO I = 1, nluns
  act_line = lineuptransitions(I)
  stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%stat_waight
  exci_energy = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
  ! internal jump up
  CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, act_popup)
  actVal = act_popup * linelist(act_line)%A_ul * exci_energy
  Zintupjump = Zintupjump + stat_weight * actVal
  Lintupjump(I) = actVal
 END DO
! the total sum 
Ztotal = Zintdownjump + Zraddeexc + Zintupjump + Zcoll
! a random number for computation, which process occurs
rand = ran2(idum)
!print*, 'do_ipackage: Ztotal = ', Ztotal
rand = rand * Ztotal
!print*, 'do_ipackage: random number: ', rand
! these variables are only to the whole line won't be too long
Z0 = Zintdownjump
Z2 = Zintupjump
Z1 = Zraddeexc
Z2 = Zcoll
!print*, 'Zintdownjump, Zintupjump, Zraddeexc: ', Zintdownjump, Zintupjump, Zraddeexc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump
! in this case a macro-atom transits into a lower state without an energy emission
IF(rand >= 0 .AND. rand < Z0) THEN
!print*, 'internal downard jump will occur'
 count_intdownjump = count_intdownjump + 1
! next transition will be an internal downward jump
 summ = 0
 DO I = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintdownjump(I)) THEN
   actual_state = linelist(linetransitions(I))%lower
   !print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
   EXIT
  END IF
  summ = summ + Lintdownjump(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! radiative deexciation
! in this case the packet is now transformed back to the r-packet
! we have to choose a new frequency, which would be calculated randomly from the 
! possible transition last_line -> some lower line
ELSE IF (rand >= Z0 .AND. rand <= Z1) THEN
!print*, 'radiative downwnard jump will occur'
count_scattering = count_scattering + 1
 summ = Z0
! next transition will be an radiative deexcitation
! print*, 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
 package(pack_index)%typ = type_rpkt
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
   IF(linelist(K)%upper == last_line) THEN
    nlns = nlns + 1
   END IF
  END IF
 END DO
 ALLOCATE(lineradtransitions(nlns), Lrad(nlns))
 J = 0
 DO K = 1, ntransitions
  IF(linelist(K)%indexe == element_index .AND. linelist(K)%indexi == ion_index) THEN
   ! transitions to a lower level
   IF(linelist(K)%upper == last_line) THEN
    J = J + 1
    lineradtransitions(J) = K
    actVal = linelist(I)%A_ul
    Lrad(J) = actVal
    Zrad = Zrad + actVal
   END IF
  END IF
 END DO
 rand = rand * Zrad
 summ = 0
 ! looking for the given line
 DO line = 1, nlns
  IF(rand >= summ .AND. rand <= summ + Lrad(line)) THEN
   ! we found the given cell now we have to compute only a new frequency
   new_freq = linelist(lineradtransitions(line))%freq
   package(pack_index)%freq_cmf = new_freq
   CALL doppler_factor(pack_index, D)
   package(pack_index)%freq_rf = package(pack_index)%freq_cmf / D
  END IF
 END DO
 active = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF (rand >= Z1 .AND. rand <= Z2) THEN
!print*, 'internal upwnard jump will occur'
 count_intupjump = count_intupjump + 1
 summ = Z1
 DO I = 1, nluns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintupjump(I)) THEN
   actual_state = linelist(lineuptransitions(I))%upper
   !print*, 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
   EXIT
  END IF
  summ = summ + Lintupjump(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional deexcitation
! in this case a macro-atom transits into a upper state without an energy emission
ELSE IF(rand >= Z2 .AND. rand <= Z3)
 summ = Z2
 DO = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lcoll(I)) THEN
   actual_state = linelist(linetransitions(I))%lower
   !print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
   ! now it will transform into a k-packet
   EXIT
  END IF
  summ = summ + Lcoll(I)
 END DO
END IF


! only one loop
!STOP 'testing the code'

DEALLOCATE(linetransitions, lineuptransitions, Lintdownjump, Lintupjump, Lraddeexc)

END DO

END SUBROUTINE do_ipackage
