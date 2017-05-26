SUBROUTINE radiative_rates(nlns, linetransitions, nluns, lineuptransitions, population, &
Lintdown, Zintdown, Lintup, Zintup, Zrad)
USE types
IMPLICIT NONE

! input variables
INTEGER                                 :: nlns, nluns
INTEGER, DIMENSION(nlns)                :: linetransitions
INTEGER, DIMENSION(nluns)               :: lineuptransitions
DOUBLE PRECISION                        :: stat_weight, exci_energy_u, exci_energy_l
DOUBLE PRECISION                        :: exci_energy
DOUBLE PRECISION                        :: population
INTEGER                                 :: element_index, ion_index
INTEGER                                 :: act_line
DOUBLE PRECISION                        :: actVal
INTEGER                                 :: I
! output variables
DOUBLE PRECISION, DIMENSION(nlns)       :: Lintdown
DOUBLE PRECISION, DIMENSION(nluns)      :: Lintup
DOUBLE PRECISION                        :: Zintdown, Zintup, Zrad

Zintdown = 0.D0
Zrad= 0.D0
Zintup = 0.D0
! we have to know which element and ion we are calculating data for
! we will use the knowledge of lines and assume that at it is
! possible at least one transition upwards or downwards and from
! the first element we get the element and the ion informations
IF(SIZE(linetransitions) /= 0) THEN
 element_index = linelist(linetransitions(1))%indexe
 ion_index = linelist(linetransitions(1))%indexi
ELSE! IF(SIZE(lineuptransitions) /= 0) THEN
 element_index = linelist(lineuptransitions(1))%indexe
 ion_index = linelist(lineuptransitions(1))%indexi
END IF 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal downward jump and radiative deexcitation
DO I = 1, nlns
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! 1.) internal downward jump
 act_line = linetransitions(I)
 ! the basic variables
 stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%stat_waight
 exci_energy_u = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%upper)%exci_energy
 exci_energy_l = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 ! internal downward jump
 ! calculation of a rate coefficient
 actVal = population * linelist(act_line)%A_ul * exci_energy_l
! print*, 'do_ipackage: stat_waight, exci_energy, population, linelist(act_line)%A_ul, actVal', &
!       stat_weight, exci_energy, population, linelist(act_line)%A_ul, actVal
 Lintdown(I) = actVal
 Zintdown = Zintdown + stat_weight * actVal
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! radiative deexcitation
 !print*, 'do_ipackage: population = ', population
 actVal = population * linelist(act_line)%A_ul * (exci_energy_u - exci_energy_l)
 Zrad = Zrad + stat_weight * actVal
! Lrad(I) = actVal
END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! internal upward jump
DO I = 1, nluns
 act_line = lineuptransitions(I)
 stat_weight = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%stat_waight
 exci_energy = elements(element_index)%ions(ion_index)%levels(linelist(act_line)%lower)%exci_energy
 ! internal jump up
 actVal = population * linelist(act_line)%A_ul * exci_energy_l
 Zintup = Zintup + stat_weight * actVal
 Lintup(I) = actVal
 !print*, 'radiative_rates: Zintup = ', Zintup, ' Lintup(I) = ', Lintup(I)
END DO


END SUBROUTINE radiative_rates
