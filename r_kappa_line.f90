SUBROUTINE r_kappa_line(pack_index, current_mgi, nextLine, nnextlines, line_dist, actirrates, tau_line)

USE types
USE constants
USE rates_r
USE dummypacket

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

DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: cur_pos, cur_dir
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: pos_min, pos_pls!, pos_lin
DOUBLE PRECISION                                :: cmf_min, cmf_pls
DOUBLE PRECISION                                :: cur_freq_rf
DOUBLE PRECISION                                :: delta
DOUBLE PRECISION                                :: deriv
DOUBLE PRECISION                                :: s_min, s_pls
DOUBLE PRECISION                                :: tau_line_2

! testing the optical depth in line calculation
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: rad_unit
DOUBLE PRECISION, DIMENSION(const_dimofspace)   :: vel_vec_1, vel_vec_2 
INTEGER                                         :: cur_dummypack, dummypack_index
                                                                        
DOUBLE PRECISION                                :: delta_r, delta_v, deriv2

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

 ! only for the testing purpose
 cur_pos = package(pack_index)%pos
 cur_dir = package(pack_index)%dir
 cur_freq_rf = package(pack_index)%freq_rf

 delta = line_dist/2.0

 s_min = line_dist - delta
 s_pls = line_dist + delta

 pos_min = cur_pos + cur_dir * s_min
 pos_pls = cur_pos + cur_dir * s_pls

 cur_dummypack = find_free_index()
 dummypack_index = cur_dummypack + SIZE(package)
 CALL copy_package(pack_index, cur_dummypack)
 ! for testing purposes with analytical homologous approximation (only the position of a packet is needed)
 dummypackage(cur_dummypack)%pos = pos_min
 CALL velo(dummypack_index, pos_min, current_mgi, vel_vec_1, 0)
 dummypackage(cur_dummypack)%pos = pos_pls
 CALL velo(dummypack_index, pos_pls, current_mgi, vel_vec_2, 0)

 rad_unit = cur_pos/norm2(cur_pos)
 
 delta_v = dot_product(vel_vec_1, rad_unit) - dot_product(vel_vec_2, rad_unit)
 delta_r = norm2(pos_min) - norm2(pos_pls)

 deriv2 = delta_v/delta_r
 write(*,*) 'r_kappa_line: vel_vec_1 = ', vel_vec_1, ' vel_vec_2 = ', vel_vec_2
 write(*,*) 'r_kappa_line: delta_v = ', delta_v, ' delta_r = ', delta_r
 write(*,*) 'r_kappa_line: deriv2 = ', deriv2

 dummypackage(cur_dummypack)%pos = pos_min
 CALL cmf_freq(dummypack_index, pos_min, cur_freq_rf, cmf_min, 0)
 dummypackage(cur_dummypack)%pos = pos_pls
 CALL cmf_freq(dummypack_index, pos_pls, cur_freq_rf, cmf_pls, 0)

 write(*,*) 'r_kappa_line: cur_dummypack = ', cur_dummypack

 deriv = abs((s_pls - s_min)/(cmf_pls - cmf_min))
 write(*,*) 'r_kappa_line: delta_s = ', s_pls - s_min, ' delta_nu = ', cmf_pls - cmf_min

 tau_line_2 = low_pop * constanta * light_speed * f_lu / (4.0 * fr_line) * corrFactor * deriv
 
 write(*,*) 'r_kappa_line: tau_line = ', actirrates%Lline(I), ' tau_line_2 = ', tau_line_2
 ! STOP 'r_kappa_line: testing'
 write(72,*) norm2(cur_pos)/R_star, actirrates%Lline(I), tau_line_2

 actirrates%nline(I) = indexline

 
 tau_line = tau_line + actirrates%Lline(I)
END DO


END SUBROUTINE
