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
INTEGER                         :: dummypackage
! calculation of optical depth
! the basic variables


constanta = (pi * e_charge**2)/( me_g * light_speed)
dummypackage = SIZE(package)
package(dummypackage) = package(pack_index)

! n_eqf_lines = 1
! DO I = nextLine + 1, ntransitions
!  IF(linelist(I)%freq == linelist(nextLine)%freq) THEN
!   n_eqf_lines = n_eqf_lines + 1
!   CYCLE
!  END IF
!  EXIT
! END DO
! nnextlines = n_eqf_lines

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
 ! write(*,*) 'r_kappa_line: el = ', exci_energy_l, ' eu = ', exci_energy_u
 Blu = light_speed**2.0 / (2.0 * h * linelist(indexline)%freq**3.0) * &
  stat_weight_u / stat_weight_l * linelist(I)%A_ul
 ! write(*,*) 'r_kappa_line: indexe = ', indexe, ' indexi = ', indexi
 ! write(*,*) 'indexl = ', linelist(indexline)%lower, 'indexu = ', linelist(indexline)%upper
 ! write(*,*) 'current_mgi = ', current_mgi
 CALL populations(indexe, indexi, linelist(indexline)%lower,&
  current_mgi, low_pop)
 CALL populations(indexe, indexi, linelist(indexline)%upper,&
  current_mgi, upp_pop)
 IF(low_pop <= 1.E-20 .OR. upp_pop == 1.E-20) THEN
  actirrates%Lline(I) = 0.E0
  actirrates%nline(I) = indexline
  CYCLE
 END IF
 ! write(*,*) 'r:r_kappa_line: indexline = ', nextLine + I - 1
 ! write(*,*) 'r_kappa_line: lower = ', linelist(indexline)%lower,&
 !  ' upper = ', linelist(indexline)%upper
 ! write(*,*) 'r_kappa_line: low_pop = ', low_pop, ' upp_pop = ', upp_pop
 corrFactor = 1.D0 - (stat_weight_l * upp_pop) / (stat_weight_u * low_pop)
 IF(corrFactor < 0.D0) THEN
  write(*,*) 'r_kappa_line: upp_pop / low_pop = ', upp_pop / low_pop
  write(*,*) 'low_pop = ', low_pop, ' upp_pop = ', upp_pop
  write(*,*) 'WARNING: correction factor 1 - (gl nu) / (gu nl) < 0'
 END IF
 !!!!!!!!!!!
 ! ROverV
 ROverV = roverw()
 ! write(*,*) 'r_kappa_line: ROverV = ', ROverV, ' low_pop = ', low_pop, &
 !  ' Blu = ', Blu, ' corrFactor = ', corrFactor, ' ldist = ', ldist / R_star, &
 !  ' freq = ', linelist(indexline)%freq
 ! write(*,*) 'r_kappa_line: rho = ', model_grid(current_mgi)%rho, ' t = ', model_grid(current_mgi)%T
 ! calculation of optical depth and rates
 ! actirrates%Lline(I) = low_pop * Blu * h * light_speed * ROverV &
 ! / (4.0 * pi) * corrFactor * ldist
 ! actirrates%Lline(I) = low_pop * pi * e_v ** 2.0  * ROverV &
 !  / (me_g * linelist(indexline)%freq) * linelist(indexline)%f_lu * corrFactor
 ! actirrates%Lline(I) = light_speed / linelist(indexline)%freq * constanta * &
 !  linelist(indexline)%f_lu * low_pop * corrFactor * ROverV
 actirrates%Lline(I) = light_speed / linelist(indexline)%freq * constanta * &
  linelist(indexline)%f_lu * low_pop * corrFactor * ROverV
 actirrates%nline(I) = indexline
 !write(*,*) 'r_kappa_line: Lline(', I, ') = ', actirrates%Lline(I)
 ! write(*,*) 'r_kappa_line: f_lu = ', linelist(indexline)%f_lu
 tau_line = tau_line + actirrates%Lline(I)
END DO

END SUBROUTINE
