SUBROUTINE r_kappa_line(pack_index, current_mgi, nextLine, nnextlines, actirrates, tau_line)

USE types
USE rates_r
IMPLICIT NONE

TYPE(rrates)                    :: actirrates
INTEGER                         :: nextLine, nnextlines, pack_index
INTEGER                         :: current_mgi
DOUBLE PRECISION                :: Blu, exci_energy_l, exci_energy_u
DOUBLE PRECISION                :: stat_weight_l, stat_weight_u
INTEGER                         :: I
INTEGER                         :: indexe, indexi, indexline
DOUBLE PRECISION                :: tau_line
INTEGER                         :: lower_level, upper_level
DOUBLE PRECISION                :: ROverV
DOUBLE PRECISION                :: low_pop, upp_pop
DOUBLE PRECISION                :: corrFactor
DOUBLE PRECISION                :: constanta
DOUBLE PRECISION                :: roverw


! the basic variables
constanta = (pi * e_charge**2)/( me_g * light_speed)

tau_line = 0.D0
DO I = 1, nnextlines
 indexline = nextLine + I - 1
 indexe = linelist(indexline)%indexe
 indexi = linelist(indexline)%indexi

 lower_level = linelist(indexline)%lower
 upper_level = linelist(indexline)%upper

 stat_weight_u = &
  elements(indexe)%ions(indexi)%levels(upper_level)%stat_waight
 stat_weight_l = &
  elements(indexe)%ions(indexi)%levels(lower_level)%stat_waight
 exci_energy_u = &
  elements(indexe)%ions(indexi)%levels(upper_level)%exci_energy
 exci_energy_l = &
  elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy

 Blu = light_speed**2.0 / (2.0 * h * linelist(indexline)%freq**3.0) * &
  stat_weight_u / stat_weight_l * linelist(I)%A_ul
 
 CALL populations(indexe, indexi, lower_level, current_mgi, low_pop)
 CALL populations(indexe, indexi, upper_level, current_mgi, upp_pop)
 
 IF(low_pop <= 1.E-20 .OR. upp_pop <= 1.E-20) THEN
  actirrates%Lline(I) = 0.E0
  actirrates%nline(I) = indexline
  CYCLE
 END IF

 corrFactor = 1.D0 - (stat_weight_l * upp_pop) / (stat_weight_u * low_pop)

 IF(corrFactor < 0.D0) THEN
  write(*,*) 'r_kappa_line: upp_pop / low_pop = ', upp_pop / low_pop
  write(*,*) 'low_pop = ', low_pop, ' upp_pop = ', upp_pop
  write(*,*) 'WARNING: correction factor 1 - (gl nu) / (gu nl) < 0'
 END IF
 !!!!!!!!!!!
 ! ROverV
 ROverV = roverw()

 actirrates%Lline(I) = light_speed / linelist(indexline)%freq * constanta * &
  linelist(indexline)%f_lu * low_pop * corrFactor * ROverV

 actirrates%nline(I) = indexline
 
 tau_line = tau_line + actirrates%Lline(I)
END DO

END SUBROUTINE
