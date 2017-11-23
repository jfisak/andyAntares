MODULE rates_k
 ! k-packages
 TYPE, PUBLIC :: krates
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcont
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcool_excit
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcool_ff
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcool_ion
  ! L(1,i): index of element, L(2,i): index of ion, L(3,i): level index
  INTEGER(KIND=2), DIMENSION(:,:), ALLOCATABLE :: Lcool_fbind
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Lcool_fbE
 END TYPE krates

 INTERFACE krates
  MODULE PROCEDURE krates_construct
 END INTERFACE

 CONTAINS
 FUNCTION krates_construct
  TYPE(krates)          :: krates_construct

  ALLOCATE(krates_construct%Lcool_excit(ntransitions))
  ALLOCATE(krates_construct%Lcool_ff(n_coll))
  ALLOCATE(krates_construct%Lcool_ion(n_phcs))
  ALLOCATE(krates_construct%Lcool_fbE(n_phcs))
 END FUNCTION krates_construct
 
 

END MODULE rates_k
