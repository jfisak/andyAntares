SUBROUTINE r_kappa_line(pack_index, current_mgi, nextLine, nnextlines, line_dist, actirrates, tau_line)

USE types
USE constants
USE rates_r
IMPLICIT NONE

INTEGER                         :: pack_index
TYPE(rrates)                    :: actirrates
INTEGER                         :: nextLine, nnextlines
DOUBLE PRECISION                :: line_dist
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
DOUBLE PRECISION                :: fr_line, f_lu


! the basic variables
constanta = (pi * e_charge**2)/( me_g * light_speed)

! write(*,*) 'r_kappa_line: nnextlines = ', nnextlines

tau_line = 0.D0
DO I = 1, nnextlines
 indexline = nextLine + I - 1
 indexe = linelist(indexline)%indexe
 indexi = linelist(indexline)%indexi

 lower_level = linelist(indexline)%lower
 upper_level = linelist(indexline)%upper

 fr_line = linelist(indexline)%freq

 stat_weight_u = &
  elements(indexe)%ions(indexi)%levels(upper_level)%stat_waight
 stat_weight_l = &
  elements(indexe)%ions(indexi)%levels(lower_level)%stat_waight
 exci_energy_u = &
  elements(indexe)%ions(indexi)%levels(upper_level)%exci_energy
 exci_energy_l = &
  elements(indexe)%ions(indexi)%levels(lower_level)%exci_energy

 Blu = light_speed**2.0 / (2.0 * h * fr_line**3.0) * &
  stat_weight_u / stat_weight_l * linelist(I)%A_ul
 
 CALL populations(indexe, indexi, lower_level, current_mgi, low_pop)
 CALL populations(indexe, indexi, upper_level, current_mgi, upp_pop)
 
 ! write(*,*) 'r_kappa_line: low_pop = ', low_pop, ' upp_pop = ', upp_pop
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
 ROverV = roverw(pack_index, line_dist, fr_line)
 f_lu = linelist(indexline)%f_lu
 actirrates%Lline(I) = light_speed / fr_line * constanta * &
  f_lu * low_pop * corrFactor * ROverV

 write(*,*) 'r_kappa_line: f_lu = ', f_lu, 'freq = ', fr_line, ' Blu = ', Blu, ' ROverV = ', ROverV

 actirrates%nline(I) = indexline
 
 tau_line = tau_line + actirrates%Lline(I)
END DO

END SUBROUTINE
