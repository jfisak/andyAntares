SUBROUTINE read_populations()

USE types
USE constants

IMPLICIT NONE

CHARACTER(80)                             :: populationfile

INTEGER                                         :: cur_element, cur_ion, cur_level
INTEGER                                         :: n_ions, n_levels

populationfile=TRIM(inputpopFile)

write(*,*) 'read_populations: populationfile = ', populationfile

! allocation of all population arrays
OPEN(UNIT=13, FILE=populationfile)
 DO cur_element = 1, n_elements
  n_ions = SIZE(elements(cur_element)%ions)
  DO cur_ion = 1, n_ions
   n_levels = SIZE(elements(cur_element)%ions(cur_ion)%levels)
   DO cur_level = 1, n_levels
    ALLOCATE(elements(cur_element)%ions(cur_ion)%levels(cur_level)%population(n_modelgrid))
    READ(13, *) elements(cur_element)%ions(cur_ion)%levels(cur_level)%population(:)
    write(*,*) 'cur_pop: cur_element = ', cur_element, ' cur_ion = ', cur_ion, ' cur_level = ', cur_level, &
     ' pop = ', elements(cur_element)%ions(cur_ion)%levels(cur_level)%population(:)
   END DO
  END DO
 END DO
CLOSE(13)

END SUBROUTINE read_populations
