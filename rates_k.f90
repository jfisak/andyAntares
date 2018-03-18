MODULE rates_k
 USE types
 ! k-packages
 TYPE, PUBLIC :: krates
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcont
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcool_excit
!  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcool_ff
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE  :: Lcool_ion
  ! L(1, i): index of element, L(2, i): index of ion, L(3, i): level index,
  ! L(4, i): cooling rate, L(5, i): starting frequency point
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Lcool_fbE
 END TYPE krates

 INTERFACE krates
  MODULE PROCEDURE krates_construct
 END INTERFACE

 CONTAINS
 FUNCTION krates_construct()
  TYPE(krates)          :: krates_construct

  ! cool_excit
  ! we expect that number of included ions does not change in the stellar wind
  ! so we do not reallocate existing array
  IF(.NOT. ALLOCATED(krates_construct%Lcool_excit)) THEN
   ALLOCATE(krates_construct%Lcool_excit(ntransitions))
  END IF

  ! cool_ff
  ! number of possible rates
  ! IF(.NOT. ALLOCATED(krates_construct%Lcool_ff)) THEN
  !  n_coll = 0
  !  DO indexe = 1, n_elements
  !   n_ions = SIZE(elements(indexe)%ions)
  !   DO indexi = 1, n_ions
  !    n_coll = n_coll + 1
  !   END DO
  !  END DO
  !  !write(*,*) 'cool_ff: n_coll = ', n_coll
  !  ! we expect that number of included ions does not change in the stellar wind
  !  ! so we do not reallocate existing array
  !  ! ALLOCATE(krates_construct%Lcool_ff(n_coll))
  ! END IF

  ! cool_ionization
  IF(.NOT. ALLOCATED(krates_construct%Lcool_ion)) THEN
   n_phcs = 0
   DO indexe = 1, n_elements
    n_ions = SIZE(elements(indexe)%ions)
    DO indexi = 1, n_ions - 1
     n_levels = SIZE(elements(indexe)%ions(indexi)%levels)
     DO indexl = 1, n_levels
      n_phcs = n_phcs + 1
     END DO ! levels
    END DO ! ions
   END DO ! elements
   ALLOCATE(krates_construct%Lcool_ion(n_phcs))
   ALLOCATE(krates_construct%Lcool_fbE(5,n_phcs))
  END IF

 END FUNCTION krates_construct

END MODULE rates_k
