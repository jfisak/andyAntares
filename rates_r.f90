MODULE rates_r

USE types
IMPLICIT NONE

!DOUBLE PRECISION, ALLOCATABLE           :: Lcont(:)
TYPE rrates
 DOUBLE PRECISION, ALLOCATABLE                   :: Lcont(:,:)

END TYPE rrates

INTERFACE rrates
  MODULE PROCEDURE rrates_construct
END INTERFACE rrates

CONTAINS
 FUNCTION rrates_construct()
  TYPE(rrates)          :: rrates_construct

  IF(.NOT. ALLOCATED(rrates_construct%Lcont)) &
   ALLOCATE(rrates_construct%Lcont(4, n_tot_cont))
!   write(*,*) 'rates_r: dim(lcont) = ', SIZE(rrates_construct%Lcont)
 END FUNCTION



END MODULE rates_r
