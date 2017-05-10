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
INTEGER                         :: nlns 
INTEGER                         :: element_index, ion_index
INTEGER                         :: I, J
INTEGER, ALLOCATABLE            :: linetransitions(:), transitions(:)
DOUBLE PRECISION, ALLOCATABLE   :: Lintdownjump(:), Lraddeexc(:)
INTEGER                         :: act_line, actVal
DOUBLE PRECISION                :: ran2, rand
! partition function
DOUBLE PRECISION                :: Z, Ztotal, Zintdownjump, Zraddeexc
! partition function for the given process
DOUBLE PRECISION                :: Z0, Z1
DOUBLE PRECISION                :: summ

! define the needed variables
last_line = package(pack_index)%last_line
last_level = linelist(last_line)%upper
element_index = linelist(last_line)%indexe
ion_index = linelist(last_line)%indexi

active = 1
! this is an initial state of the macro-atom
actual_state = last_level
! we will run this loop until the macro atom is deactivated
DO WHILE (active == 1)
 ! we have to find all possible downward upward transitions
 ! firstly we calculate number of these possible transitions
 DO I = 1, ntransitions
  ! we are interested only in the transitions for the given atom
  IF(linelist(I)%indexe == element_index .AND. linelist(I)%indexi == ion_index) THEN
   IF(linelist(I)%upper == actual_state) THEN
    nlns = nlns + 1
   END IF
  END IF
 END DO
 
 ALLOCATE(linetransitions(nlns))
 ! we will save these possible transitions into an array
 J = 0
 DO I = 1, ntransitions
  IF(linelist(I)%indexe == element_index .AND. linelist(I)%indexi == ion_index) THEN
   IF(linelist(I)%upper == actual_state) THEN
    J = J + 1
    linetransitions(J) = I
   END IF
  END IF
 END DO
 
 ! total rate of an internal downward jump within the current
 ! ionization state
 Zintdownjump = 0.D0
 ALLOCATE(Lintdownjump(nlns))
 Zraddeexc = 0.D0
 ALLOCATE(Lraddeexc(nlns))
 DO I = 1, nlns
  act_line = linetransitions(I)
  ! internal jump down
  actVal = linelist(act_line)%A_ul
  Zintdownjump = Zintdownjump + actVal
  Lintdownjump(I) = actVal
  ! radiative deexcitation
  actVal = linelist(act_line)%A_ul * &
   (elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy - &
   elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy)
  Zraddeexc = Zraddeexc + actVal
  Lraddeexc(I) = actVal
 END DO
! the total statistical sum 
Ztotal = Zintdownjump + Zraddeexc
! a random number for computation, which process occurs
rand = ran2(idum) * Z
! these variables are only to the whole line won'ŧ be too long
Z0 = Zintdownjump
Z1 = Zraddeexc
IF(rand < Z0) THEN
! next transition will be an internal downward jump
 summ = 0
 DO I = 1, nlns
  ! we will find the given state
  IF(rand >= summ .AND. rand < summ + Lintdownjump(I)) actual_state = linelist(linetransitions(I))%lower
  summ = summ + Lintdownjump(I)
 END DO
ELSE IF (rand <= Z0 .AND. rand <= Z1) THEN
 summ = Z1
 DO I = 1, nlns
! next transition will be an radiative deexcitation
  IF(rand >= summ .AND. rand < summ + Lraddeexc(I)) THEN
   summ = summ + Lraddeexc(I)
   package(pack_index)%typ = type_rpkt
   CALL emit_rpackage(pack_index)
   active = 0
  END IF
 END DO
END IF




END DO

END SUBROUTINE do_ipackage
