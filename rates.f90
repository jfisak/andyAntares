MODULE rates

! pointers
! i-packages
! internal radiative downward jump
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_dorad => NULL()
! internal upward jump
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_uprad => NULL()
! radiative deexcitation
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_rad => NULL()
! internal collisional downward jump
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_docoll => NULL()
! internal collisional upward jump
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_upcoll => NULL()
! internal upward jump
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_up => NULL()
! internal downward jump
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_do => NULL()
! radiative recombination
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_recrad => NULL()
! internal radiative jump to the lower ionization state
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_recrad => NULL()
! collisional recombination
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_reccol => NULL()
! internal collisional jump to the lower ionization state
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_reccol => NULL()

! k-packages
DOUBLE PRECISION, DIMENSION(:), POINTER  ::  Lcont => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_excit => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_ff => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_ion => NULL()
! L(1,i): index of element, L(2,i): index of ion, L(3,i): level index
INTEGER(KIND=2), DIMENSION(:,:), POINTER :: Lcool_fbind => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_fbE => NULL()

END MODULE rates
