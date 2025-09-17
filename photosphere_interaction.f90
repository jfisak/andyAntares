SUBROUTINE photosphere_interaction(pack_index)

! this sbr decides what happens with packet if
! it flies back to the photosphere and it is not
! destroyed
USE types
USE constants
IMPLICIT NONE

INTEGER                                         :: pack_index, n_pack
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: directionn, direction
DOUBLE PRECISION                                :: doppler_D, L_star
DOUBLE PRECISION                                :: sint, cost, sinp, cosp
INTEGER                                         :: ind_cell_numb
DOUBLE PRECISION, DIMENSION(const_dimofspace)                  :: corner, width, pos, upcorner
DOUBLE PRECISION                                :: freq
n_pack = SIZE(package) - 1
L_star = 4.D0*const_pi*(R_star)**2*const_stefbolz*T_eff**4
SELECT CASE(abs_surface)
! creation of a new packet
CASE(1)
 CALL random_unitvector1(direction, sint, cost, sinp, cosp)
 ! write(*,*) 'init_photsphere: R_star = ', R_star
 package(pack_index)%pos = R_star * direction

 ! Then give it a random direction outward from the photosphere
 CALL random_unitvector2(directionn) !random_unitvector(direction)
 direction(ind_x)=directionn(ind_z)*sint*cosp+directionn(ind_x)*cost*cosp-directionn(ind_y)*sinp
 direction(ind_y)=directionn(ind_z)*sint*sinp+directionn(ind_x)*cost*sinp+directionn(ind_y)*cosp
 direction(ind_z)=directionn(ind_z)*cost-directionn(ind_x)*sint
 package(pack_index)%dir = direction
 CALL find_dyn_cell1(package(pack_index)%pos,ind_cell_numb)
 IF(ind_cell_numb > SIZE(dyn_cell)) THEN
  write(*,*) 'init_photsphere: wrong cell number'
  CALL abort()
 END IF
 ! write(*,*) 'photosphere_interaction: ind_cell_numb = ', ind_cell_numb
 package(pack_index)%cell_numb = ind_cell_numb
 corner = dyn_cell(ind_cell_numb)%corner
 upcorner = dyn_cell(ind_cell_numb)%upcorner
 width = dyn_cell(ind_cell_numb)%width
 pos = package(pack_index)%pos
 
 CALL freq_from_planck(freq, T_eff)   ! here the frequency is sampled from the Planck law
 package(pack_index)%freq_rf = freq
 package(pack_index)%e_rf = L_star/n_pack  
 CALL doppler_factor(pack_index, doppler_D)
 package(pack_index)%freq_cmf = package(pack_index)%freq_rf * doppler_D 
 package(pack_index)%e_cmf    = package(pack_index)%e_rf * doppler_D  
 package(pack_index)%last_line = no_line
 package(pack_index)%delta_s = 0.D0

 package(pack_index)%n_interactions = 0
 package(pack_index)%next_cross = NONE


IF(pos(ind_x) < corner(ind_x) .OR. pos(ind_x) > upcorner(ind_x) .OR. &
 pos(ind_y) < corner(ind_y) .OR. pos(ind_y) > upcorner(ind_y) .OR. &
 pos(ind_z) < corner(ind_z) .OR. pos(ind_z) > upcorner(ind_z)) THEN
  write(*,*) 'find_dist: pack_index = ', pack_index
  write(*,*) 'find_dist: pos/corner = ', pos(:)/corner(:)!, ' corner = ', corner / R_inf
  write(*,*) 'find_dist: pos/R_inf = ', pos/R_inf, ' corner/R_inf = ', corner/R_inf
  write(*,*) 'find_dist: cell_numb = ', ind_cell_numb, ' neighbors = ', dyn_cell(ind_cell_numb)%neighbor
  STOP 
END IF

CASE DEFAULT
 write(99,*) 'photosphere_interaction: the choice'
 write(99,*) abs_surface, ' is not possible'
 STOP
END SELECT

END SUBROUTINE photosphere_interaction
