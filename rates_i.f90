MODULE rates_i

! pointers
! i-packages
! pre-calculated values of integrals
! array of temperatures
DOUBLE PRECISION, ALLOCATABLE                   :: i_temps(:)
TYPE i_ints
 INTEGER                                        :: indexe, indexi, indexl
 DOUBLE PRECISION, ALLOCATABLE                  :: gammaijk(:), alphaijk(:)
 ! DOUBLE PRECISION, ALLOCATABLE                  :: etaijk(:)
END TYPE
TYPE(i_ints), ALLOCATABLE                       :: iints(:)

! internal radiative downward jump
TYPE, PUBLIC :: irates
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_dorad 
 ! internal upward jump
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_uprad 
 ! radiative deexcitation
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_rad 
 ! internal collisional downward jump
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_docoll 
 ! internal collisional upward jump
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_upcoll 
 ! internal upward jump
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_up 
 ! internal downward jump
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_do 
 ! radiative recombination
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_recrad 
 ! internal radiative jump to the lower ionization state
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_recrad 
 ! collisional recombination
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_reccol 
 ! internal collisional jump to the lower ionization state
 DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lma_int_reccol 
 ! destructor
 !SUBROUTINE, PASS :: rates_destruct
END TYPE irates

INTERFACE irates
 module procedure rates_construct
END INTERFACE

CONTAINS
FUNCTION rates_construct(nlns, nluns, nlio)
 TYPE(irates)                   :: rates_construct
 ! INTEGER                        :: nldo, nlup
 INTEGER                        :: nlio

 ALLOCATE(rates_construct%Lma_int_dorad(nlns))
 ALLOCATE(rates_construct%Lma_int_docoll(nlns))
 ALLOCATE(rates_construct%Lma_rad(nlns))
 ALLOCATE(rates_construct%Lma_int_uprad(nluns))
 ALLOCATE(rates_construct%Lma_int_upcoll(nluns))
 ALLOCATE(rates_construct%Lma_int_do(nlns))
 ALLOCATE(rates_construct%Lma_int_up(nluns))
 ALLOCATE(rates_construct%Lma_recrad(nlio))
 ALLOCATE(rates_construct%Lma_int_recrad(nlio))
 ALLOCATE(rates_construct%Lma_reccol(nlio))
 ALLOCATE(rates_construct%Lma_int_reccol(nlio))
END FUNCTION rates_construct



END MODULE rates_i
