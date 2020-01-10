DOUBLE PRECISION FUNCTION roverv()

USE types
IMPLICIT NONE

DOUBLE PRECISION                               :: R_pos


IF(velapprox == 0) THEN
 ROverV = R_inf / V_inf
 ! actirrates%Lline(I) = low_pop * Blu * h * light_speed * ROverV &
 ! / (4.0 * pi) * corrFactor * ldist
 ! write(*,*) 'r_kappa_line: actirrates%Lline(I) = ', actirrates%Lline(I)
ELSE IF(velapprox == 2) THEN
 ! according to (10) in Abbot & Lucy (1985)
 ! r
 R_pos = norm2(package(dummypackage)%pos)
 ! ||v||
 V_pos = V_inf * (1.0 - R_star / R_pos ) ** beta
 ! v = (v_x, v_y, v_z)
 V_pos_vec = V_pos * package(dummypackage)%pos / norm2(package(dummypackage)%pos)
 costheta = dot_product(package(dummypackage)%dir, V_pos_vec) / norm2(V_pos_vec)
 ROverV = (V_inf - V_0) / (R_inf - R_star) + 1 / R_pos * (1 - costheta**2.0) *&
  (V_0 * R_inf - V_inf * R_star) / (R_inf - R_star)
ELSE IF(velapprox == 1) THEN
 ! package(dummypackage) = package(pack_index)
 ! CALL emit_rpackage(dummypackage)
 CALL boundary3(pack_index, cell_dist, next_cell)
 ! according to (10) in Abbot & Lucy (1985)
 ! r
 R_pos = norm2(package(dummypackage)%pos)
 ! ||v||
 V_pos = V_inf * (1.0 - R_star / R_pos ) ** beta
 ! v = (v_x, v_y, v_z)
 V_pos_vec = V_pos * package(pack_index)%pos / norm2(package(pack_index)%pos)
 ! \mu
 costheta = dot_product(package(pack_index)%dir, V_pos_vec) / norm2(V_pos_vec)
 ! dv/dr
 dV_pos = beta * R_star * V_inf / R_pos**2 * (1.0 - R_star / R_pos)**(beta-1)
 ROverV = 1.0 / (costheta**2.0 * dV_pos + (1.0 - costheta**2.0)* V_pos / R_pos)
 ! actirrates%Lline(I) = low_pop * Blu * h * light_speed * &
 !  ROverV / (4.0 * pi) * corrFactor 
END IF

return ROverV

END FUNCTION
