SUBROUTINE saha_factor(indexe, indexi, indexl, temp, sfactor)

USE types

IMPLICIT NONE

DOUBLE PRECISION                        :: temp
INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: dimind
! INTEGER, ALLOCATABLE                    :: indexl0(:)
INTEGER, DIMENSION(1)                   :: indexl0
DOUBLE PRECISION                        :: popj, popjp1
DOUBLE PRECISION                        :: sfactor
DOUBLE PRECISION                        :: sahaconst, gijk, g0, eijk, e0

!CALL populations(indexe, indexi - 1, indexl, cur_mgi, popj)
!! dimind = SIZE(MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy))
!! ALLOCATE(indexl0(dimind))
indexl0(:) = MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
!CALL populations(indexe, indexi, indexl0, cur_mgi, popjp1)
!
!sfactor = popj / (popjp1 * electron_density)
!
!write(*,*) 'saha_factor: indexi = ', indexi, ' indexl = ', indexl
!write(*,*) 'saha_factor: popj = ', popj, ' popjp1 = ', popjp1, &
! ' electron_density = ', electron_density
!write(*,*) 'saha_factor: sfactor = ', sfactor

sahaconst = 5.D-1 * (h**2 / (2.0 * pi * me_g * BOLK))**(3.0/2.0)
gijk = elements(indexe)%ions(indexi - 1)%levels(indexl)%stat_waight
g0 = elements(indexe)%ions(indexi)%levels(indexl0(1))%stat_waight
eijk = elements(indexe)%ions(indexi - 1)%levels(indexl)%exci_energy
e0 = elements(indexe)%ions(indexi)%levels(indexl0(1))%exci_energy

sfactor = gijk / g0 * sahaconst / temp**(3.0/2.0) * exp((e0 - eijk) / (BOLK * temp))

! write(*,*) 'saha_factor: gijk = ', gijk, ' g0 = ', g0, ' temp = ', temp,&
!  ' e0 = ', e0, ' eijk = ', eijk
! write(*,*) 'saha_factor: sfactor = ', sfactor

END SUBROUTINE saha_factor
