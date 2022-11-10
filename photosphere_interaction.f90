SUBROUTINE photosphere_interaction(pack_index)

! this sbr decides what happens with packet if
! it flies back to the photosphere and it is not
! destroyed
USE types
IMPLICIT NONE

INTEGER                                         :: pack_index, n_pack
DOUBLE PRECISION, DIMENSION(3)                  :: directionn, direction
DOUBLE PRECISION                                :: D, L_star
DOUBLE PRECISION                                :: sint, cost, sinp, cosp
INTEGER                                         :: ind_cell_numb
DOUBLE PRECISION, DIMENSION(3)                  :: corner, width, pos
DOUBLE PRECISION                                :: freq
n_pack = SIZE(package) - 1
L_star = 4.D0*pi*(R_star)**2*sigma*T_eff**4
SELECT CASE(abs_surface)
! creation of a new packet
CASE(1)
 CALL random_unitvector1(direction, sint, cost, sinp, cosp)
 ! write(*,*) 'init_photsphere: R_star = ', R_star
 package(pack_index)%pos = R_star * direction

 ! Then give it a random direction outward from the photosphere
 CALL random_unitvector2(directionn) !random_unitvector(direction)
 direction(1)=directionn(3)*sint*cosp+directionn(1)*cost*cosp-directionn(2)*sinp
 direction(2)=directionn(3)*sint*sinp+directionn(1)*cost*sinp+directionn(2)*cosp
 direction(3)=directionn(3)*cost-directionn(1)*sint
 package(pack_index)%dir = direction
 CALL find_dyn_cell1(package(pack_index)%pos,ind_cell_numb)
 IF(ind_cell_numb > SIZE(dyn_cell)) THEN
  write(*,*) 'init_photsphere: wrong cell number'
  CALL abort()
 END IF
 ! write(*,*) 'photosphere_interaction: ind_cell_numb = ', ind_cell_numb
 package(pack_index)%cell_numb = ind_cell_numb
 corner = dyn_cell(ind_cell_numb)%corner
 width = dyn_cell(ind_cell_numb)%width
 pos = package(pack_index)%pos
 
 CALL freq_from_planck(freq, T_eff)   ! here the frequency is sampled from the Planck law
 package(pack_index)%freq_rf = freq
 package(pack_index)%e_rf = L_star/n_pack  
 CALL doppler_factor(pack_index, D)
 package(pack_index)%freq_cmf = package(pack_index)%freq_rf * D 
 package(pack_index)%e_cmf    = package(pack_index)%e_rf * D  
 package(pack_index)%last_line = no_line
 package(pack_index)%delta_s = 0.D0

 package(pack_index)%n_interactions = 0
 package(pack_index)%next_cross = NONE


IF(pos(1) < corner(1) .OR. pos(1) > corner(1) + width(1) .OR. &
 pos(2) < corner(2) .OR. pos(2) > corner(2) + width(2) .OR. &
 pos(3) < corner(3) .OR. pos(3) > corner(3) + width(3) ) THEN
  write(*,*) 'find_dist: pack_index = ', pack_index
  write(*,*) 'find_dist: pos/corner = ', pos(:)/corner(:)!, ' corner = ', corner / R_inf
  write(*,*) 'find_dist: pos/R_inf = ', pos/R_inf, ' corner/R_inf = ', corner/R_inf,&
   ' width/R_inf = ', width/R_inf
  write(*,*) 'find_dist: cell_numb = ', ind_cell_numb, ' neighbors = ', dyn_cell(ind_cell_numb)%neighbor
  STOP 
END IF

CASE DEFAULT
 write(99,*) 'photosphere_interaction: the choice'
 write(99,*) abs_surface, ' is not possible'
 STOP
END SELECT

END SUBROUTINE photosphere_interaction
