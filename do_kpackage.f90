SUBROUTINE do_kpackage(pack_index)
USE types
IMPLICIT NONE

INTEGER                         :: pack_index
! indexes
INTEGER                         :: I
! cooling rates
DOUBLE PRECISION                :: Zexcit, Ztot
DOUBLE PRECISION                :: Z0
DOUBLE PRECISION, DIMENSION(ntransitions) &
                                :: Lcool_excit, Lexcit
DOUBLE PRECISION                :: summ
DOUBLE PRECISION                :: rand, ran2
!package(pack_index)%typ = type_rpkt

! calculating of cooling rates
! collision excitation rate
CALL cool_excit(1, pack_index,  Zexcit, Lexcit)


! ještě to bude dobré seřadit od nejpravděpodobnější po nejméně pravděpodobný proces

! now we have to decide which cooling process will occure
rand = ran2(idum)
Ztot = Zexcit
rand = rand * Ztot
Z0 = Zexcit
!print*, 'do_kpackage: Ztot = ', Ztot, ' rand = ', rand
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! collisional excitation
IF(rand >= 0.D0 .AND. rand <= Z0) THEN
 summ = 0.D0
 DO I = 1, ntransitions
  !print*, 'Lexcit(I) = ', Lexcit(I)
  IF(rand >= summ .AND. rand < summ + Lexcit(I)) THEN
   !print*, 'do_kpackage: packet = ', pack_index, ' a collisional excitation occures'
   !print*, 'collisional deexcitation: I = ', I, ' upper level = ', linelist(I)%upper
   package(pack_index)%last_line = I
   package(pack_index)%typ = type_ipkt
   EXIT
  END IF
  summ = summ + Lexcit(I)
 END DO
!STOP 'do_kpackage: testing'
END IF



END SUBROUTINE do_kpackage
