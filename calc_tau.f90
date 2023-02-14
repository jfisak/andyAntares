! a subroutine to calculate optical depth from the particular point and direction
SUBROUTINE calc_tau(init_pos, direction, frequency)

USE types
USE rates_r
IMPLICIT NONE


DOUBLE PRECISION, DIMENSION(3)                  :: init_pos, direction
DOUBLE PRECISION                                :: frequency
LOGICAL                                         :: active
INTEGER                                         :: cur_packet, cur_approx
DOUBLE PRECISION                                :: cur_s

INTEGER                                         :: next_cell, init_line, next_line, n_lines, cur_cell
INTEGER                                         :: cur_mgi, get_package_model_index
DOUBLE PRECISION                                :: bound_dist, line_dist, dopf
DOUBLE PRECISION                                :: kappa_cont, kappa, tau_line
LOGICAL                                         :: in_cell, cell_change

TYPE(rrates)                                    :: actirrates

DOUBLE PRECISION                                :: cur_tau_line, cur_tau_cont, cur_r

INTEGER                                         :: n_thomson

active = .true.
cur_packet = 1
cur_approx = 1

cur_s = 0.D0
cur_tau_line = 0.D0
cur_tau_cont = 0.D0

! number of thomson scattering
! total number of continuum opacity sources
IF(n_ff /= 0) THEN
! nothing to be done in this fork
ELSE
  n_ff = 1
  n_thomson = 1
  n_tot_cont = n_thomson + n_photcrossect + n_ff
END IF
actirrates = rrates()


CALL find_dyn_cell1(init_pos, cur_cell)
package(cur_packet)%cell_numb = cur_cell
package(cur_packet)%pos = init_pos
package(cur_packet)%dir = direction
package(cur_packet)%last_line = no_line
package(cur_packet)%next_cross = NONE
package(cur_packet)%freq_rf = frequency
CALL doppler_factor(cur_packet, dopf)
package(cur_packet)%freq_cmf = frequency * dopf

init_line = no_line

OPEN(148, file="tau.dat")

! while the packet is active
DO WHILE(active)
 ! distance to the cell boundary
 CALL boundary3(cur_packet, bound_dist, next_cell)
 ! find next line
 ! write(*,*) 'calc_tau: next_line_bluered'
 CALL next_line_bluered(cur_approx, cur_packet, bound_dist, init_line, next_line, n_lines)
 ! write(*,*) 'calc_tau: next_line = ', next_line
 ALLOCATE(actirrates%Lline(n_lines), actirrates%nline(n_lines))
 ! find a distance to the next resonance point
 CALL resonance_distance2(cur_packet, next_line, bound_dist, in_cell, line_dist)
 
 ! continuum kappa
 CALL r_kappa_cont(cur_packet, kappa_cont, actirrates)
 CALL doppler_factor(cur_packet, dopf)
 ! line tau
 cur_mgi = get_package_model_index(cur_packet)
 CALL r_kappa_line(cur_packet, cur_mgi, next_line, n_lines, line_dist, actirrates, tau_line)
 ! write(*,*) 'calc_tau: tau_line = ', tau_line, ' line_dist = ', line_dist, dopf * frequency / linelist(next_line)%freq


 if(line_dist < bound_dist) then
  cell_change = .FALSE.
  CALL move_package(cur_packet, line_dist, next_cell, cell_change)
  
  cur_r = norm2(package(cur_packet)%pos)
  cur_s = cur_s + line_dist
  cur_tau_cont = cur_tau_cont + dopf * kappa_cont * line_dist
  write(148, *) cur_r, cur_s, cur_tau_line, cur_tau_cont

  cur_tau_line = cur_tau_line + tau_line
  write(148, *) cur_r, cur_s, cur_tau_line, cur_tau_cont

 else if (line_dist > bound_dist .and. bound_dist > 0.e0) then
  cell_change = .TRUE.
  CALL move_package(cur_packet, bound_dist, next_cell, cell_change)

  cur_r = norm2(package(cur_packet)%pos)
  cur_s = cur_s + bound_dist
  cur_tau_cont = cur_tau_cont + dopf * kappa_cont * bound_dist
  write(148, *) cur_r, cur_s, cur_tau_line, cur_tau_cont
 else if (bound_dist < 0.e0) then
  CALL change_cell(cur_packet, next_cell)
 end if

 if(cur_r > R_inf .or. cur_r < R_star) EXIT 

 init_line = next_line

 DEALLOCATE(actirrates%Lline, actirrates%nline)
 
END DO

CLOSE(148)

END SUBROUTINE calc_tau
