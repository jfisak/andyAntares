! calculation of the opacity in the lines
! which lines? the first line which Sobolev point is the closest to the path of the packet
!
! INPUT: pack_index(INT): index of the packet
!        current_mgi(INT): index of the modGrid point
!        nextLine(INT): index of the first line (in the array)
!        nnextlines(INT): number of lines with the same frequency as the nextLine
! OUTPUT: line_dist(DBLE): a distance to the Sobolev point
!        actirrates(rrates): rates saved for specific transitions
!        tau_line(DBLE): optical depth in the chosen lines
!
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
DOUBLE PRECISION                                :: tau_line_2, tau_line_3

INTEGER, PARAMETER                              :: max_n_of_velopackets = 20000

! testing the optical depth in line calculation
INTEGER                                         :: cur_dummypack, dummypack_index

! the basic variables
constanta = (const_pi * const_e**2)/( const_me_g * const_c)

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
 f_lu = linelist(indexline)%f_lu
 IF(velApprox == 0 .or. velApprox == 1) THEN
  ROverV = roverw(pack_index, line_dist, fr_line)
  actirrates%Lline(I) = const_c / fr_line * constanta * &
   f_lu * low_pop * corrFactor * ROverV
  actirrates%nline(I) = indexline
 ELSE IF(velApprox == 3) THEN


  cur_pos = package(pack_index)%pos
  cur_dir = package(pack_index)%dir
  cur_freq_rf = package(pack_index)%freq_rf
 
  delta = basic_cell_width(1)/100.0
 
  s_min = line_dist - delta
  s_pls = line_dist + delta
 
  pos_min = cur_pos + cur_dir * s_min
  pos_pls = cur_pos + cur_dir * s_pls
 
  cur_dummypack = find_free_index()
  dummypack_index = cur_dummypack + SIZE(package)
  CALL copy_package(pack_index, cur_dummypack)
  ! for testing purposes with analytical homologous approximation (only the position of a packet is needed)
  ! move packet to the pos_min
  ! CALL teleport_dummypacket(cur_dummypack, pos_min)
  ! CALL velo(dummypack_index, vel_vec_1, 0)
  ! ! move packet to the pos_pls
  ! CALL teleport_dummypacket(cur_dummypack, pos_pls)
  ! CALL velo(dummypack_index, vel_vec_2, 0)
 
  ! rad_unit1 = pos_pls/norm2(pos_pls)
  ! rad_unit2 = pos_min/norm2(pos_min)
  ! 
  ! delta_v = dot_product(vel_vec_1, rad_unit1) - dot_product(vel_vec_2, rad_unit2)
  ! delta_r = norm2(pos_min) - norm2(pos_pls)
 
  ! IF(delta_v == 0.D0) THEN
  !  deriv2 = 0.D0
  ! ELSE
  !  deriv2 = delta_r/delta_v
  ! END IF
  ! write(*,*) 'r_kappa_line: vel_vec_1 = ', vel_vec_1, ' vel_vec_2 = ', vel_vec_2
  ! write(*,*) 'r_kappa_line: delta_v = ', delta_v, ' delta_r = ', delta_r
  ! write(*,*) 'r_kappa_line: deriv2 = ', deriv2
  ! write(73,*) norm2(cur_pos)/R_star, deriv2
 
  ! move packet to the pos_min
  CALL teleport_dummypacket(cur_dummypack, pos_min)
  CALL cmf_freq(dummypack_index, cur_freq_rf, cmf_min)
  ! move packet to the pos_pls
  CALL teleport_dummypacket(cur_dummypack, pos_pls)
  CALL cmf_freq(dummypack_index, cur_freq_rf, cmf_pls)
 
  ! write(*,*) 'r_kappa_line: cur_dummypack = ', cur_dummypack
 
  deriv = abs((s_pls - s_min)/(cmf_pls - cmf_min))
  ! write(*,*) 'r_kappa_line: delta_s = ', s_pls - s_min, ' delta_nu = ', cmf_pls - cmf_min
 
  tau_line = low_pop * constanta * f_lu * corrFactor * deriv
  ! write(*,*) 'r_kappa_line: tau_line = ', tau_line
  ! tau_line_2 = const_c / fr_line * constanta * f_lu * low_pop * corrFactor * deriv2
  actirrates%Lline(I) = tau_line
  actirrates%nline(I) = indexline

  ! only for the testing purpose
  ! costheta = dot_product(package(pack_index)%dir, V_pos_vec) / V_pos
  ! ROverV = 1.0 / (costheta**2.0 * dV_pos + (1.0 - costheta**2.0)* V_pos / R_pos)
  ! ROverV = roverw(pack_index, line_dist, fr_line)
  tau_line_3 = const_c / fr_line * constanta * &
   f_lu * low_pop * corrFactor * R_inf/V_inf
  ! actirrates%Lline(I) = tau_line
  IF(pack_index < max_n_of_velopackets) THEN
   write(72,*) norm2(cur_pos)/R_star, tau_line_3, tau_line_2, deriv
  END IF
 
  CALL deactivate_dummy_packet(cur_dummypack)
 
 END IF

 tau_line = tau_line + actirrates%Lline(I)
END DO


END SUBROUTINE
