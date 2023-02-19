! a subroutine to calculate optical depth from the particular point and direction
SUBROUTINE calc_tau(init_pos, end_pos, frequency)

USE types
USE rates_r
IMPLICIT NONE


DOUBLE PRECISION, DIMENSION(3)                  :: init_pos, end_pos
DOUBLE PRECISION, DIMENSION(3)                  :: direction, cur_pos, cross_pos
DOUBLE PRECISION                                :: frequency
INTEGER                                         :: next_cross
LOGICAL                                         :: active, snapped
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

DOUBLE PRECISION, PARAMETER                     :: largeNumber = 1.D99
DOUBLE PRECISION                                :: txm, txp, typ, tym, tzm, tzp

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
 ELSE IF(cur_r > R_inf .or. cur_r < R_star) THEN
  cur_pos = package(cur_packet)%pos
  txp = 0.D0
  txm = 0.D0
  typ = 0.D0
  tym = 0.D0
  tzp = 0.D0
  tzm = 0.D0

  txm = (cur_pos(1) - xmin)/direction(1)
  txp = (cur_pos(1) - xmax)/direction(1)
  tym = (cur_pos(2) - ymin)/direction(2)
  typ = (cur_pos(2) - ymax)/direction(2)
  tzm = (cur_pos(3) - zmin)/direction(3)
  tzp = (cur_pos(3) - zmax)/direction(3)

  if(txm < bound_dist .and. txm > 0.D0) then
   next_cross = negx
   bound_dist = txm
  end if
  if(txp < bound_dist .and. txp > 0.D0) then
   next_cross = posx
   bound_dist = txp
  end if
  if(tym < bound_dist .and. tym > 0.D0) then
   next_cross = negy
   bound_dist = tym
  end if
  if(typ < bound_dist .and. typ > 0.D0) then
   next_cross = posy
   bound_dist = typ
  end if
  if(tzm < bound_dist .and. txm > 0.D0) then
   next_cross = negz
   bound_dist = tzm
  end if
  if(tzp < bound_dist .and. txp > 0.D0) then
   next_cross = posz
   bound_dist = tzp
  end if

  cross_pos = cur_pos + direction * bound_dist
  if(next_cross == txm .or. next_cross == txp) then
   ! 
   if(cross_pos(2) > ymin .and. cross_pos(2) < ymax .and. &
    cross_pos(3) > zmin .and. cross_pos(3) < zmax) then
     snapped = .true.
   end if
  end if
  
  if(next_cross == tym .or. next_cross == typ) then
   if(cross_pos(1) > xmin .and. cross_pos(1) < xmax .and. &
    cross_pos(3) > zmin .and. cross_pos(3) < zmax) then
     snapped = .true.
   end if
   
  end if
  if(next_cross == tzm .or. next_cross == tzp) then
   if( cross_pos(2) > ymin .and. cross_pos(2) < ymax .and. &
    cross_pos(1) > xmin .and. cross_pos(1) < xmax) then
     snapped = .true.
   end if
   
  end if

  if(snapped) then
   package(1)%pos = package(1)%pos + bound_dist * package(1)%dir 
  else
   exit
  end if
  
 END IF

 init_line = next_line

 DEALLOCATE(actirrates%Lline, actirrates%nline)
 
END DO

CLOSE(148)

END SUBROUTINE calc_tau
