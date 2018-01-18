SUBROUTINE do_kpackage(pack_index)
USE types
USE rates_k
IMPLICIT NONE

INTEGER                         :: pack_index
! indexes
INTEGER                         :: I
! cooling rates
DOUBLE PRECISION                :: Zexcit, Ztot, Zff, Zion, Zfb
DOUBLE PRECISION                :: Z0, Z1, Z2, Z3
DOUBLE PRECISION                :: summ
DOUBLE PRECISION                :: rand
!package(pack_index)%typ = type_rpkt
! total number of possible cooling processes
INTEGER                         :: n_cool_tot
! new frequency
DOUBLE PRECISION                :: new_freq
! choosing the given process
INTEGER                         :: act_proc, J
TYPE(krates)                    :: actikrates
REAL(8)                         :: random
! indexes for ion levels
INTEGER                         :: indexe, indexi, indexl
INTEGER                         :: actIndex
INTEGER                         :: n_ions, n_levels
INTEGER                         :: nline

actikrates = krates()

! calculating of cooling rates
! collision excitation rate
CALL cool_excit(1, pack_index,  Zexcit, actikrates)
CALL cool_ff(pack_index, Zff, actikrates)
CALL cool_ionization(1, pack_index, Zion, actikrates)
CALL cool_fb(pack_index, Zfb, actikrates)

! now we allocate an array which will include all possible transitions
! n_cool_tot = SIZE(Lcool_excit) + SIZE(Lcool_ff)
! ALLOCATE(cool_rates(n_cool_tot))
! adding processes into fields
! 1.) cooling excitations

! now we have to decide which cooling process will occure
rand = random()
Z0 = Zexcit
Z1 = Z0 + Zff
Z2 = Z1 + Zion
Z3 = Z2 + Zfb
! total rate
Ztot = Zexcit + Zff + Zion + Zfb
rand = rand * Ztot
! write(*,*) 'do_kpackage: Zexcit = ', Zexcit, ' Zff = ', Zff, ' Zion = ', Zion, ' Zfb = ', Zfb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional excitation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(rand >= 0.D0 .AND. rand <= Z0) THEN
 summ = 0.D0
 DO I = 1, ntransitions
  !print*, 'Lcool_excit(I) = ', Lcool_excit(I)
  IF(rand >= summ .AND. rand < summ + actikrates%Lcool_excit(I)) THEN
   !print*, 'do_kpackage: packet = ', pack_index, ' a collisional excitation occures'
   !print*, 'collisional deexcitation: I = ', I, ' upper level = ', linelist(I)%upper
   package(pack_index)%last_line = I
   package(pack_index)%typ = type_ipkt
   ! write(*,*) 'do_kpackage: package = ', pack_index, ' collisional excitation process...'
   EXIT
  END IF
  summ = summ + actikrates%Lcool_excit(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FREE-FREE processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF(rand >= Z0 .AND. rand <= Z1) THEN
 summ = Z0
 ! free-free process
 ! package changes to r-packet
 package(pack_index)%typ = type_rpkt
 package(pack_index)%last_line = no_line
 CALL k_freq_ff(pack_index, new_freq)
 write(*,*) 'do_kpackage: package = ', pack_index, ' free-free process...'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ionization
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF(rand >= Z1 .AND. rand <= Z2) THEN
 summ = Z1
 package(pack_index)%typ = type_ipkt
 actIndex = 0
 DO indexe = 1, n_elements
  n_ions = SIZE(elements(indexe)%ions)
  DO indexi = 1, n_ions - 1
   n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
   DO indexl = 1, n_levels
    IF(rand > summ .AND. rand <= summ + actikrates%Lcool_ion(actIndex + 1)) THEN
     ! we have to find any corresponding line for the given ion including the given
     ! line in the lower or the upper level(, which must be saved as well)
     DO nline = 1, ntransitions
      IF(indexe == linelist(actIndex)%indexe .AND. &
         indexi == linelist(actIndex)%indexi) THEN
       IF(indexl == linelist(actIndex)%upper) THEN
        package(pack_index)%last_line = nline
        isUpperTransition = .TRUE.
        EXIT
       ELSE IF(indexl == linelist(actIndex)%lower) THEN
        package(pack_index)%last_line = nline
        isUpperTransition = .FALSE.
        EXIT
       ! test for level number
       END IF
      ! test for indexe and indexi
      END IF
     ! loop over lines
     END DO
    END IF
   ! loop over levels
   END DO
  ! loop over ions
  END DO
 ! loop over elements
 END DO
 write(*,*) 'do_kpackage: package = ', pack_index, ' ionization process...'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! recombination
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF(rand > Z2 .AND. rand <= Z3) THEN
 summ = Z2
 !write(*,*) 'do_kpackage: rand = ', rand, 'Z2 = ', Z2, ' Z3 = ', Z3
 DO J = 1, SIZE(actikrates%Lcool_fbE(:))
  !write(*,*) 'do_kpackage: rand = ', rand, 'summ = ', summ, ' summ + Lcool_fbE = ', &
  ! summ + Lcool_fbE(J)
  IF(rand > summ .AND. rand <= summ + actikrates%Lcool_fbE(J)) THEN
   act_proc = J
   write(*,*) 'do_kpackage: act_proc = ', act_proc
   EXIT
  END IF
  summ = summ + actikrates%Lcool_fbE(J)
 END DO
 ! package changes to r-packet
 package(pack_index)%typ = type_rpkt
 package(pack_index)%last_line = no_line
 write(*,*) 'do_kpackage: package = ', pack_index, ' recombination process...'
 CALL k_freq_fb(pack_index, act_proc, new_freq)



!STOP 'do_kpackage: testing'
END IF



END SUBROUTINE do_kpackage
