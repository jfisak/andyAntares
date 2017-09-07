SUBROUTINE do_kpackage(pack_index)
USE types
USE rates
IMPLICIT NONE

INTEGER                         :: pack_index
! indexes
INTEGER                         :: I
! cooling rates
DOUBLE PRECISION                :: Zexcit, Ztot, Zff
DOUBLE PRECISION                :: Z0, Z1
DOUBLE PRECISION                :: summ
DOUBLE PRECISION                :: rand, ran2
!package(pack_index)%typ = type_rpkt
! total number of possible cooling processes
INTEGER                         :: n_cool_tot
! new frequency
DOUBLE PRECISION                :: new_freq


! calculating of cooling rates
! collision excitation rate
CALL cool_excit(1, pack_index,  Zexcit)
CALL cool_ff(pack_index, Zff)

! now we allocate an array which will include all possible transitions
!n_cool_tot = SIZE(Lcool_excit) + SIZE(Lcool_ff)
!ALLOCATE(cool_rates(n_cool_tot))
! adding processes into fields
! 1.) cooling excitationskk

! now we have to decide which cooling process will occure
rand = ran2(idum)
Z0 = Zexcit
Z1 = Z0 + Zff
! total rate
Ztot = Zexcit + Zff
rand = rand * Ztot
print*, 'do_kpackage: Zexcit = ', Zexcit, ' Zff = ', Zff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional excitation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(rand >= 0.D0 .AND. rand <= Z0) THEN
 summ = 0.D0
 DO I = 1, ntransitions
  !print*, 'Lcool_excit(I) = ', Lcool_excit(I)
  IF(rand >= summ .AND. rand < summ + Lcool_excit(I)) THEN
   !print*, 'do_kpackage: packet = ', pack_index, ' a collisional excitation occures'
   !print*, 'collisional deexcitation: I = ', I, ' upper level = ', linelist(I)%upper
   package(pack_index)%last_line = I
   package(pack_index)%typ = type_ipkt
   EXIT
  END IF
  summ = summ + Lcool_excit(I)
 END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FREE-FREE processes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF(rand >= Z0 .AND. rand <= Z1) THEN
 summ = Z1
 ! free-free process
 ! package changes to r-packet
 package(pack_index)%typ = type_rpkt
 package(pack_index)%last_line = no_line
 CALL k_freq_ff(pack_index, new_freq)
 write(*,*) 'do_kpackage: package = ', pack_index, ' free-free process...'



!STOP 'do_kpackage: testing'
END IF



END SUBROUTINE do_kpackage
