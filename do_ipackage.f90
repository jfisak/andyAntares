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
INTEGER                         :: I, J, K
INTEGER, ALLOCATABLE            :: linetransitions(:), transitions(:), lineuptransitions(:)
DOUBLE PRECISION, ALLOCATABLE   :: Lintdownjump(:), Lraddeexc(:), Lintupjump(:)
INTEGER                         :: act_line, actVal
DOUBLE PRECISION                :: ran2, rand
! partition function
DOUBLE PRECISION                :: Z, Ztotal, Zintdownjump, Zraddeexc, Zintupjump
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1, Z2
DOUBLE PRECISION                :: summ, stat_weight, exci_energy
! populations
DOUBLE PRECISION                :: act_popup, act_popdown
INTEGER                         :: get_package_model_index, current_mgi

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
 
 ! total rate of an internal downward jump within the current
 ! ionization state
 Zintdownjump = 0.D0
 ALLOCATE(Lintdownjump(nlns))
 Zintupjump = 0.D0
 ALLOCATE(Lintupjump(nluns))
 Zraddeexc = 0.D0
 ALLOCATE(Lraddeexc(nlns))
 DO I = 1, nlns
  act_line = linetransitions(I)
  ! internal jump down
  ! calculation of number density of the given ion
  CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, act_popdown)
  ! statistical weight
  stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
  exci_energy = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
  ! calculation of a rate coefficient
  !print*, 'do_ipackage: stat_waight, exci_energy, act_popdown', stat_weight, exci_energy, act_popdown
  actVal = act_popdown * linelist(act_line)%A_ul * exci_energy/e_v
  Zintdownjump = Zintdownjump + stat_weight * actVal
  Lintdownjump(I) = actVal
  ! radiative deexcitation
  CALL populations(element_index, ion_index, linelist(act_line)%upper, current_mgi, act_popup)
  !print*, 'do_ipackage: act_popup = ', act_popup
  actVal = act_popup * linelist(act_line)%A_ul * &
   (elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy - &
   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy)
  Zraddeexc = Zraddeexc + stat_weight * actVal
  Lraddeexc(I) = actVal
 END DO
!print*, 'do_ipackage: Zintdownjump = ', Zintdownjump
 DO I = 1, nluns
  act_line = lineuptransitions(I)
  ! internal jump up
  actVal = linelist(act_line)%A_ul
  Zintupjump = Zintupjump + actVal
  Lintupjump(I) = actVal
 END DO
! the total statistical sum 
Ztotal = Zintdownjump + Zraddeexc + Zintupjump
! a random number for computation, which process occurs
rand = ran2(idum)
!print*, 'do_ipackage: Ztotal = ', Ztotal
rand = rand * Ztotal
!print*, 'do_ipackage: random number: ', rand
! these variables are only to the whole line won'ŧ be too long
Z0 = Zintdownjump
Z2 = Zintupjump
Z1 = Zraddeexc
IF(rand < Z0) THEN
! next transition will be an internal downward jump
 summ = 0
 DO I = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintdownjump(I)) THEN
   actual_state = linelist(linetransitions(I))%lower
   print*, 'do_ipackage: packet: ', pack_index, ' internal downward jump...'
   EXIT
  END IF
  summ = summ + Lintdownjump(I)
 END DO
ELSE IF (rand <= Z0 .AND. rand <= Z1) THEN
 summ = Z0
 DO I = 1, nlns
! next transition will be an radiative deexcitation
  IF(rand >= summ .AND. rand < summ + Lraddeexc(I)) THEN
 print*, 'do_ipackage: packet: ', pack_index, ' radiative deexcitation...'
   summ = summ + Lraddeexc(I)
   package(pack_index)%typ = type_rpkt
   CALL emit_rpackage(pack_index)
   active = 0
  END IF
 END DO
ELSE IF (rand <= Z1 .AND. rand <= Z2) THEN
 summ = Z1
 DO I = 1, nluns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintdownjump(I)) THEN
   actual_state = linelist(linetransitions(I))%upper
   print*, 'do_ipackage: packet: ', pack_index, ' internal upward jump...'
   EXIT
  END IF
  summ = summ + Lintdownjump(I)
 END DO
END IF

! only one loop
EXIT

DEALLOCATE(linetransitions, lineuptransitions, Lintdownjump, Lintupjump, Lraddeexc)

END DO

END SUBROUTINE do_ipackage
