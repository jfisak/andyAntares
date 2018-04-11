SUBROUTINE r_kappa_line(pack_index, current_mgi, nextLine, nnextlines, ldist, Lline)

USE types
IMPLICIT NONE

INTEGER                         :: nextLine, nnextlines, pack_index
INTEGER                         :: current_mgi
DOUBLE PRECISION                :: ldist
DOUBLE PRECISION                :: Blu, exci_energy_l, exci_energy_u
DOUBLE PRECISION                :: stat_weight_l, stat_weight_u
INTEGER                         :: I
DOUBLE PRECISION, DIMENSION(nnextlines)         :: Lline
INTEGER                         :: indexe, indexi, indexline
DOUBLE PRECISION                :: R_res, V_res
DOUBLE PRECISION, DIMENSION(3)  :: V_res_vec
DOUBLE PRECISION                :: tau_line, dV_res
DOUBLE PRECISION                :: costheta
INTEGER                         :: lower_level
DOUBLE PRECISION                :: ROverV
DOUBLE PRECISION                :: low_pop, upp_pop
DOUBLE PRECISION                :: corrFactor
 
! calculation of optical depth
! the basic variables
tau_line = 0.D0
DO I = 1, nnextlines
 indexline = nextLine + I - 1
 indexe = linelist(indexline)%indexe
 indexi = linelist(indexline)%indexi
 lower_level = linelist(indexline)%lower
 stat_weight_u = &
  elements(indexe)%ions(indexi)%levels(linelist(indexline)%upper)%stat_waight
 stat_weight_l = &
  elements(indexe)%ions(indexi)%levels(linelist(indexline)%lower)%stat_waight
 exci_energy_u = &
  elements(indexe)%ions(indexi)%levels(linelist(indexline)%upper)%exci_energy
 exci_energy_l = &
  elements(indexe)%ions(indexi)%levels(linelist(indexline)%lower)%exci_energy
 Blu = light_speed**2.0 / (2.0 * h * linelist(indexline)%freq**3.0) * &
  stat_weight_u / stat_weight_l * linelist(I)%A_ul
 CALL populations(indexe, indexi, linelist(nextLine + I - 1)%lower,&
  current_mgi, low_pop)
 CALL populations(indexe, indexi, linelist(nextLine + I - 1)%upper,&
  current_mgi, upp_pop)
 write(*,*) 'r:r_kappa_line: indexline = ', nextLine + I - 1
 write(*,*) 'r_kappa_line: lower = ', linelist(nextLine + I - 1)%lower,&
  ' upper = ', linelist(nextLine + I - 1)%upper
 write(*,*) 'r_kappa_line: low_pop = ', low_pop
 IF(velapprox == 0) THEN
  R_res = R_inf
  V_res = V_inf
  corrFactor = 1.D0 - (stat_weight_l * upp_pop) / (stat_weight_u * low_pop)
  IF(corrFactor < 0.D0) write(*,*) 'WARNING: correction factor 1 - (gl nu) / (gu nl) < 0'
  Lline(I) = low_pop * Blu * h * light_speed * (R_res / V_res) &
  / (4.0 * pi) * corrFactor
  
  write(*,*) 'resonance_distance: low_pop = ', low_pop, ' Blu = ', Blu, &
   'stat_weight_l = ', stat_weight_l, ' stat_weight_u = ', stat_weight_u
  write(*,*) 'resonance_distance: upp_pop = ', upp_pop,&
   ' stat_weight_l = ', stat_weight_l, ' stat_weight_u = ', stat_weight_u
  write(*,*) 'resonance_distance: 1 - (gl nu) / (gu nl) = ',&
   (1.D0 - (stat_weight_l * upp_pop) / (stat_weight_u * low_pop))
  write(*,*) 'resonance_distance: Lline(I) = ', Lline(I)
 ELSE IF(velapprox == 1) THEN
  R_res = norm2(package(pack_index)%pos + package(pack_index)%dir * ldist)
  V_res = V_inf * (1.0 - R_star / R_res ) ** beta
  costheta = dot_product(package(pack_index)%dir, V_res_vec) / V_res
  dV_res = beta * R_star * V_res / R_res * (1.0 - R_star / R_res)**(-1)
  ROverV = 1.0 / (costheta**2.0 * dV_res + (1.0 - costheta**2.0)* V_res / R_res)
 END IF
 tau_line = tau_line + Lline(I)
END DO

END SUBROUTINE
