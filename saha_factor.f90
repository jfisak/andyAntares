SUBROUTINE saha_factor(indexe, indexi, indexl, cur_mgi, electron_density, sfactor)

USE types

IMPLICIT NONE

INTEGER                                 :: indexe, indexi, indexl
INTEGER                                 :: cur_mgi, dimind
! INTEGER, ALLOCATABLE                    :: indexl0(:)
INTEGER, DIMENSION(1)                   :: indexl0
DOUBLE PRECISION                        :: popj, popjp1
DOUBLE PRECISION                        :: electron_density
DOUBLE PRECISION                        :: sfactor

CALL populations(indexe, indexi - 1, indexl, cur_mgi, popj)
! dimind = SIZE(MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy))
! ALLOCATE(indexl0(dimind))
indexl0(:) = MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
! STOP 'saha_factor: testing'
CALL populations(indexe, indexi, indexl0, cur_mgi, popjp1)

sfactor = popj / (popjp1 * electron_density)

 write(*,*) 'saha_factor: indexi = ', indexi, ' indexl = ', indexl
 write(*,*) 'saha_factor: popj = ', popj, ' popjp1 = ', popjp1, &
  ' electron_density = ', electron_density
 write(*,*) 'saha_factor: sfactor = ', sfactor

END SUBROUTINE saha_factor
