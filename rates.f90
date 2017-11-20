MODULE rates

! k-packages
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcont => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_excit => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_ff => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_ion => NULL()
! L(1,i): index of element, L(2,i): index of ion, L(3,i): level index
INTEGER(KIND=2), DIMENSION(:,:), POINTER :: Lcool_fbind => NULL()
DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lcool_fbE => NULL()
! pointers
! i-packages
! internal radiative downward jump
TYPE, PUBLIC :: irates
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_dorad 
 ! internal upward jump
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_uprad 
 ! radiative deexcitation
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_rad 
 ! internal collisional downward jump
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_docoll 
 ! internal collisional upward jump
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_upcoll 
 ! internal upward jump
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_up 
 ! internal downward jump
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_do 
 ! radiative recombination
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_recrad 
 ! internal radiative jump to the lower ionization state
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_recrad 
 ! collisional recombination
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_reccol 
 ! internal collisional jump to the lower ionization state
 DOUBLE PRECISION, DIMENSION(:), POINTER  :: Lma_int_reccol 
 ! destructor
 !SUBROUTINE, PASS :: rates_destruct
END TYPE irates

TYPE(irates)      :: actirates

INTERFACE irates
 module procedure rates_construct
END INTERFACE

CONTAINS
FUNCTION rates_construct(nlns, nluns, nlio)
 TYPE(irates)                   :: rates_construct
 INTEGER                        :: nldo, nlup, nlio
 INTEGER                         :: OMP_GET_THREAD_NUM, my_rank

 my_rank = OMP_GET_THREAD_NUM()
 write(*,*) 'construction: my_rank = ', my_rank, ' nlns = ', nlns, ' nluns = ', nluns, ' nlio = ', nlio
 ALLOCATE(rates_construct%Lma_int_dorad(nlns))
 ALLOCATE(rates_construct%Lma_int_docoll(nlns))
 ALLOCATE(rates_construct%Lma_int_uprad(nluns))
 ALLOCATE(rates_construct%Lma_int_upcoll(nluns))
 ALLOCATE(rates_construct%Lma_int_do(nlns))
 ALLOCATE(rates_construct%Lma_int_up(nluns))
 ALLOCATE(rates_construct%Lma_recrad(nlio))
 ALLOCATE(rates_construct%Lma_int_recrad(nlio))
 ALLOCATE(rates_construct%Lma_reccol(nlio))
 ALLOCATE(rates_construct%Lma_int_reccol(nlio))
END FUNCTION rates_construct


!FUNCTION rates_destruct(nlns, nluns, nlio)
! TYPE(irates)                   :: rates_destruct
! INTEGER                        :: nldo, nlup, nlio
! DEALLOCATE(rates_destruct%Lma_int_dorad(nlns))
! DEALLOCATE(rates_destruct%Lma_int_docoll(nlns))
! DEALLOCATE(rates_destruct%Lma_int_uprad(nluns))
! DEALLOCATE(rates_destruct%Lma_int_upcoll(nluns))
! DEALLOCATE(rates_destruct%Lma_int_do(nlns))
! DEALLOCATE(rates_destruct%Lma_int_up(nluns))
! DEALLOCATE(rates_destruct%Lma_recrad(nlio))
! DEALLOCATE(rates_destruct%Lma_int_recrad(nlio))
! DEALLOCATE(rates_destruct%Lma_reccol(nlio))
! DEALLOCATE(rates_destruct%Lma_int_reccol(nlio))
! DEALLOCATE(rates_destruct%Lma_int_dorad(nlns))
!END FUNCTION rates_destruct
 


END MODULE rates
