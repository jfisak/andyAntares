SUBROUTINE velo(pack_index, pack_pos, vel_vec)

USE types

IMPLICIT NONE    

INTEGER                           :: pack_index
DOUBLE PRECISION                  :: vel_radial, vec_length
DOUBLE PRECISION, DIMENSION(3)    :: vel_vec, pack_pos

SELECT CASE(velApprox)
! homologous expansion
CASE(0)
 vel_radial = V_inf/R_inf * vec_length(pack_pos)
 vel_vec = pack_pos/vec_length(pack_pos) * vel_radial
! the beta velocity law
CASE(1)
 vel_radial = V_inf * (1.D0 - R_star / norm2(pack_pos))**beta
 vel_vec = pack_pos/vec_length(pack_pos) * vel_radial
CASE(2)
 vel_radial = (V_inf - V_0)/(R_inf - R_star) * &
  vec_length(pack_pos) + &
  (V_0 * R_inf - V_inf * R_star) / (R_inf - R_star)
 vel_vec = pack_pos/vec_length(pack_pos) * vel_radial
! #03
! velocity field given by model in discrete points
CASE(3)
 CALL vel_discrete_points(pack_index, vel_vec)
CASE DEFAULT
 write(*,*) 'velo: velApprox = ', velApprox
 write(*,*) 'this velocity structure is not known'
 CALL abort()
END SELECT
! check if the packet is located inside the model grid
! IF(vec_length(package(pack_index)%pos) > R_inf .OR. &
!  norm2(package(pack_index)%pos) < R_star) THEN
!  vel_vec = (/ 0.0, 0.0, 0.0/)
! END IF
 IF(norm2(vel_vec) > light_speed) THEN
  write(*,*) 'velo: V_inf = ', V_inf, ' R_star = ', R_star, ' beta = ', beta,&
   ' ||package(pack_index)%pos|| = ', norm2(pack_pos)
  write(*,*) 'velo: pack_index = ', pack_index, ' position = ', vec_length(pack_pos)/R_inf,&
  norm2(pack_pos)/R_star
  write(*,*) 'velocity is larger than the speed of light'
  CALL abort()
 END IF
END SUBROUTINE velo
