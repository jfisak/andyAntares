SUBROUTINE saha_factor(indexe, indexi, indexl, cur_mgi, electron_density, sfactor)

USE types

IMPLICIT NONE

INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: cur_mgi, dimind
! INTEGER, ALLOCATABLE                    :: indexl0(:)
INTEGER, DIMENSION(1)                   :: indexl0
DOUBLE PRECISION                        :: popj, popjp1
DOUBLE PRECISION                        :: electron_density, temp
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
temp = model_grid(cur_mgi)%t

sfactor = gijk / g0 * sahaconst / temp**(3.0/2.0) * exp((e0 - eijk) / (BOLK * temp))

!write(*,*) 'saha_factor: sfactor = ', sfactor

END SUBROUTINE saha_factor
