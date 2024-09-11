! a subroutine to calculate optical depth between two points r1 and r2
SUBROUTINE calc_tau(init_pos, end_pos, frequency)

USE types
USE rates_r
USE constants
IMPLICIT NONE


DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: init_pos, end_pos
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: direction, cur_pos
DOUBLE PRECISION                                :: frequency
LOGICAL                                         :: active
INTEGER                                         :: cur_packet, cur_approx
DOUBLE PRECISION                                :: cur_s

INTEGER                                         :: next_cell, init_line, next_line, n_lines, cur_cell
INTEGER                                         :: cur_mgi, get_package_model_index
DOUBLE PRECISION                                :: bound_dist, line_dist, dopf, end_dist
DOUBLE PRECISION                                :: kappa_cont, tau_line
LOGICAL                                         :: in_cell, cell_change

TYPE(rrates)                                    :: actirrates

DOUBLE PRECISION                                :: cur_tau_line, cur_tau_cont, cur_r

INTEGER                                         :: n_thomson


DOUBLE PRECISION                                :: distance

active = .true.
cur_packet = 1
cur_approx = 1
cur_s = 0.D0
cur_tau_line = 0.D0
cur_tau_cont = 0.D0


distance = norm2(end_pos - init_pos)
if(distance /= 0.0) then
 direction = (end_pos - init_pos)/distance
else
 write(99,*) 'calc_tau: distance is equal to zero'
 return
end if

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
 cur_pos = package(cur_packet)%pos
 cur_r = norm2(cur_pos)
 write(*,*) 'calc_tau: r = ', cur_r/R_star, cur_r/R_inf
 ! the current point is located inside the propGrid
 IF(cur_r < R_inf .or. cur_r > R_star) THEN
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

  ! distance to the end point
  end_dist = norm2(cur_pos - end_pos)
 
  write(*,*) 'calc_tau: ld = ', line_dist, ' bd = ', bound_dist, ' ed = ', end_dist
  if(line_dist < end_dist .and. line_dist < bound_dist) then
   cell_change = .false.
   CALL move_package(cur_packet, end_dist, next_cell, cell_change)
   active = .false.
  end if
  if(line_dist < bound_dist .and. active) then
   cell_change = .FALSE.
   CALL move_package(cur_packet, line_dist, next_cell, cell_change)
   
   cur_r = norm2(package(cur_packet)%pos)
   cur_s = cur_s + line_dist
   cur_tau_cont = cur_tau_cont + dopf * kappa_cont * line_dist
   write(148, *) cur_r, cur_s, cur_tau_line, cur_tau_cont
 
   cur_tau_line = cur_tau_line + tau_line
   write(148, *) cur_r, cur_s, cur_tau_line, cur_tau_cont
 
  else if (line_dist > bound_dist .and. bound_dist > 0.e0 .and. active) then
   cell_change = .TRUE.
   CALL move_package(cur_packet, bound_dist, next_cell, cell_change)
 
   cur_r = norm2(package(cur_packet)%pos)
   cur_s = cur_s + bound_dist
   cur_tau_cont = cur_tau_cont + dopf * kappa_cont * bound_dist
   write(148, *) cur_r, cur_s, cur_tau_line, cur_tau_cont
  else if (bound_dist < 0.e0 .and. active) then
   CALL change_cell(cur_packet, next_cell)
  end if
 ELSE IF(cur_pos(ind_x) < xmin .or. cur_pos(ind_x) > xmax .or. &
  cur_pos(ind_y) < ymin .or. cur_pos(ind_y) > ymax .or. &
  cur_pos(ind_z) < zmin .or. cur_pos(ind_z) > zmax) THEN
  write(*,*) 'calc_tau: calling propgrid_dist'
  CALL propgrid_dist(package(cur_packet), bound_dist)
 ELSE IF(cur_r < R_star) THEN
  active = .false.
 ELSE IF(cur_r >= R_inf) THEN
  cell_change = .false.
  CALL move_package(cur_packet, end_dist, next_cell, cell_change)
 END IF
 init_line = next_line

 DEALLOCATE(actirrates%Lline, actirrates%nline)
 
END DO

CLOSE(148)

END SUBROUTINE calc_tau
