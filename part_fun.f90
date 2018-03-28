SUBROUTINE part_fun(indexe, indexi, temp, U)

! Calculate partition function of element indexe in th egiven ionization
! stage indexi at the given temperature temp
USE types

IMPLICIT NONE    

INTEGER             :: indexe, indexi, indexl, nlevels
DOUBLE PRECISION    :: U, temp, g_level, e_level, e_gl
INTEGER, DIMENSION(1)   :: indexl0, indexl1

  
!  print*, 'partition func. called for:', indexe, indexi, temp

IF(.NOT. ALLOCATED(elements(indexe)%ions(indexi)%levels)) THEN
 write(*,*) 'part_fun: indexe = ', indexe, ' indexi = ', indexi
 write(*,*) 'levels are not allocated'
 STOP
END IF

indexl0(:) = MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
indexl1(:) = MAXLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
U = elements(indexe)%ions(indexi)%levels(indexl0(1))%stat_waight
e_gl = elements(indexe)%ions(indexi)%levels(indexl0(1))%exci_energy
!  print*, '  Part.func. initialisation:', U, e_gl
! IF(indexe == 3 .AND. indexi == 3) write(*,*) 'part_fun: indexl0(1) = ', indexl0(1), ' e_gl = ', e_gl / e_V

! Number of levels for the given element indexe in ionisation stage indexi
nlevels = SIZE(elements(indexe)%ions(indexi)%levels)

DO indexl = 1, nlevels
 IF(indexl == indexl0(1)) CYCLE
 ! Statistical weight of th egrpund level
 g_level = elements(indexe)%ions(indexi)%levels(indexl)%stat_waight   
 ! Excitation energy of the excited level
 e_level = elements(indexe)%ions(indexi)%levels(indexl)%exci_energy
 ! Partition function
 IF(temp == 0 ) STOP 'part_fun: temperature = 0...'
 U = U + g_level * EXP(-(e_level - e_gl) / BOLK / temp)  
 IF(indexe == 3 .AND. indexi == 2) write(*,*) '   part.func. calculation:', indexl, g_level, U,&
  e_level/e_V, e_gl/e_V
END DO

!  print*, '   part.func. calculation done:', U

END SUBROUTINE part_fun
