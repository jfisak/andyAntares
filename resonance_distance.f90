SUBROUTINE resonance_distance(pack_index, f_line, ldist)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: pack_index
DOUBLE PRECISION                :: f_line
! output variables
DOUBLE PRECISION                :: ldist
INTEGER                         :: approx = 0

! Calculate distance the photon needs to travel to come to
! resonance with the next line. This assumes homologous
! expansion i.e. velocity is proportional to r. Projected
! gradient of the projected velocity in a direction of the
! photon propagation is more complicated in case of no
! homologous expansion
SELECT CASE(approx)
! homologous approximation
CASE(0)
! This is case when we assume that he have only hydrogen 
ldist = light_speed * (R_inf/V_inf) * ((package(pack_index)%freq_cmf - f_line)/package(pack_index)%freq_rf)
!print*, 'resonance_distance: ldist = ', ldist
CASE DEFAULT
 STOP 'resonance_distance: this approximation is not known'
END SELECT

END SUBROUTINE resonance_distance
