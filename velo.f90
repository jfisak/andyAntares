
SUBROUTINE velo(pack_index,vel_vec,approx)

USE types

IMPLICIT NONE    

INTEGER                           :: approx
INTEGER                           :: pack_index
DOUBLE PRECISION                  :: vel_radial, vec_length
DOUBLE PRECISION, DIMENSION(3)    :: vel_vec, pack_position
! Petr Kurfurst's disk model variables
INTEGER                           :: pack_mi, get_package_model_index
DOUBLE PRECISION, DIMENSION(3)    :: vel_rad, vel_ang
DOUBLE PRECISION                  :: vel_rad_norm, vel_ang_norm

INTEGER                            :: cur_mgi

pack_position = package(pack_index)%pos
SELECT CASE(velApprox)
! homologous expansion
CASE(0)
 vel_radial = V_inf/R_inf * vec_length(pack_position)
 vel_vec = pack_position/vec_length(pack_position) * vel_radial
 write(*,*) 'velo: V_inf = ', V_inf, ' R_inf = ', R_inf
! the beta velocity law
CASE(1)
 vel_radial = V_inf * (1.D0 - R_star / norm2(pack_position))**beta
 vel_vec = pack_position/vec_length(pack_position) * vel_radial
CASE(2)
 cur_mgi = get_package_model_index(pack_index)
 vel_radial = model_grid(cur_mgi)%vel
 vel_vec = pack_position/vec_length(pack_position) * vel_radial
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
! IF(vec_length(pack_position) > R_inf .OR. &
!  norm2(pack_position) < R_star) THEN
!  vel_vec = (/ 0.0, 0.0, 0.0/)
! END IF
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! petr kurfurst's disk model
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF ((model_type .EQ. 2) .AND. (inputmodel .EQ. 1)) THEN
  ! we have to know the velocity of matter in the given point
  pack_mi = get_package_model_index(pack_index)
  vel_rad_norm = model_grid(pack_mi)%vel
  vel_ang_norm = model_grid(pack_mi)%velang
  ! now we compute given vectors
  vel_rad = (/ vel_rad_norm * pack_position(1) / vec_length(pack_position(1)), &
             vel_rad_norm * pack_position(2) / vec_length(pack_position(1)), 0.D0 /)
  vel_ang = (/ - vel_ang_norm * pack_position(2) / vec_length(pack_position(1)), &
             vel_ang_norm * pack_position(1) / vec_length(pack_position(1)), 0.D0 /)
  ! and finally the velocity vector
  vel_vec = vel_rad + vel_ang
 END IF
 ! write(*,*) 'velo: norm2(vel_vec) = ', norm2(vel_vec)
 IF(norm2(vel_vec) > light_speed) THEN
  write(*,*) 'velo: vel_vec/c = ', norm2(vel_vec)/light_speed, ' Rinf/c = ', V_inf/light_speed
  write(*,*) 'velo: V_inf = ', V_inf, ' R_star = ', R_star, ' beta = ', beta,&
   ' ||pack_position|| = ', norm2(pack_position)
  write(*,*) 'velo: pack_index = ', pack_index, ' position = ', vec_length(pack_position)/R_inf,&
  norm2(pack_position)/R_star
  write(*,*) 'velocity is larger than the speed of light'
  CALL abort()
 END IF

write(30,*) norm2(pack_position), norm2(vel_vec)

END SUBROUTINE velo
