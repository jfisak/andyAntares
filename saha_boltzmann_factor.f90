! Calculate the ratio sb_factor between population of any two exication levels 
! indexi and (indexi+1)
!
! INPUT: indexe(INT): element index
!        indexi(INT): ion index
!        temp(DBL): temperature
! OUTPUT: sb_factor(DBL): saha-boltzmann factor
!         too_large(LOG): the number is too large
!
SUBROUTINE saha_boltzmann_factor(indexe, indexi, temp, sb_factor, too_large)

USE types  
USE constants

IMPLICIT NONE

INTEGER                              :: indexe, indexi
DOUBLE PRECISION                     :: temp, sb_factor, part_U1, part_U2
DOUBLE PRECISION, PARAMETER          :: large_number = 1.D150
LOGICAL                              :: too_large

!  print*, 'Saha fact. called for:',indexe, indexi, temp

CALL part_fun(indexe, indexi, temp, part_U1)
CALL part_fun(indexe, indexi+1, temp, part_U2)
write(*,*) '  Part.func:', part_U1, part_U2, indexi
! write(*,*) 'saha_boltzmann_factor: ion_pot = ', elements(indexe)%ions(indexi)%ion_potential
sb_factor = part_U1 / part_U2 * saha_const * temp**(-3.D0/2.D0) * &
 EXP( elements(indexe)%ions(indexi)%ion_potential / (const_kB * temp) )
IF(sb_factor > large_number) THEN
 too_large = .TRUE.
ELSE
 too_large = .FALSE.
END IF

END SUBROUTINE saha_boltzmann_factor

