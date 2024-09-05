SUBROUTINE velo(pack_index, vel_vec, approx)

USE types
USE constants
USE dummypacket

IMPLICIT NONE    

INTEGER                           :: approx
INTEGER                           :: pack_index
DOUBLE PRECISION                  :: vel_radial, vec_length
DOUBLE PRECISION, DIMENSION(3)    :: vel_vec, pack_position
! Petr Kurfurst's disk model variables
INTEGER                           :: pack_mi, get_package_model_index
DOUBLE PRECISION, DIMENSION(3)    :: vel_rad, vel_ang
DOUBLE PRECISION                  :: vel_rad_norm, vel_ang_norm
DOUBLE PRECISION                  :: r_pos

LOGICAL, PARAMETER                :: velocityTesting = .true.
INTEGER, PARAMETER                :: max_n_of_velopackets = 200

INTEGER                            :: cur_mgi, cur_dummy_index


! dosti nelogické využívání dvakrát té stejné proměnné...
IF(pack_index <= SIZE(package)) THEN
 pack_position = package(pack_index)%pos
ELSE IF(pack_index > SIZE(package)) THEN
 cur_dummy_index = pack_index - SIZE(package)
 pack_position = dummypackage(cur_dummy_index)%pos
END IF

SELECT CASE(velApprox)
! #00
! homologous expansion
CASE(0)
 vel_radial = V_inf/R_inf * norm2(pack_position)
 vel_vec = pack_position/norm2(pack_position) * vel_radial
! #01
! the beta velocity law
CASE(1)
 r_pos = norm2(pack_position)
 vel_radial = V_inf * (1.D0 - R_star / norm2(pack_position))**beta
 ! write(*,*) 'velo: pack_position = ', pack_position
 if(r_pos <= R_star .or. r_pos > R_inf) vel_radial = 0.D0
 if(isnan(vel_radial)) STOP 'velo: vel_radial = NaN'
 vel_vec = pack_position/NORM2(pack_position) * vel_radial
! #02
CASE(2)
 CALL vel_discrete_points(pack_index, vel_vec)
! #03
! velocity field given by model in discrete points
CASE(3)
 ! IF(dyn_cell == 0) THEN
 CALL vel_discrete_points(pack_index, vel_vec)
 ! write(*,*) 'velo: vel_vec = ', vel_vec
 ! ELSE IF(dyn_cell /= 0) THEN
  ! calc velocity based on precalculated interpolated velocity profiles
 ! END IF
CASE DEFAULT
 write(*,*) 'velo: velApprox = ', velApprox
 write(*,*) 'this velocity structure is not known'
 CALL abort()
END SELECT

 IF(norm2(pack_position) < R_star) THEN
  vel_vec = (/ 0.D0, 0.D0, 0.D0 /)
 ELSE IF(norm2(pack_position) > R_inf) THEN
  vel_vec = V_inf * pack_position/norm2(pack_position)
 END IF
IF(velocityTesting) THEN
 IF(pack_index <= max_n_of_velopackets) THEN
  write(34,*) norm2(pack_position)/R_star, norm2(vel_vec)
 END IF
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


END SUBROUTINE velo




 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! petr kurfurst's disk model
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  IF ((model_type .EQ. 2) .AND. (inputmodel .EQ. 1)) THEN
!   ! we have to know the velocity of matter in the given point
!   pack_mi = get_package_model_index(pack_index)
!   vel_rad_norm = model_grid(pack_mi)%vel
!   vel_ang_norm = model_grid(pack_mi)%velang
!   ! now we compute given vectors
!   vel_rad = (/ vel_rad_norm * pack_position(1) / vec_length(pack_position(1)), &
!              vel_rad_norm * pack_position(2) / vec_length(pack_position(1)), 0.D0 /)
!   vel_ang = (/ - vel_ang_norm * pack_position(2) / vec_length(pack_position(1)), &
!              vel_ang_norm * pack_position(1) / vec_length(pack_position(1)), 0.D0 /)
!   ! and finally the velocity vector
!   vel_vec = vel_rad + vel_ang
!  END IF
! IF ((model_type .EQ. 2) .AND. (inputmodel .EQ. 1)) THEN
!  ! we have to know the velocity of matter in the given point
!  pack_mi = get_package_model_index(pack_index)
!  vel_rad_norm = model_grid(pack_mi)%vel
!  vel_ang_norm = model_grid(pack_mi)%velang
!  ! now we compute given vectors
!  vel_rad = (/ vel_rad_norm * pack_position(1) / vec_length(pack_position(1)), &
!             vel_rad_norm * pack_position(2) / vec_length(pack_position(1)), 0.D0 /)
!  vel_ang = (/ - vel_ang_norm * pack_position(2) / vec_length(pack_position(1)), &
!             vel_ang_norm * pack_position(1) / vec_length(pack_position(1)), 0.D0 /)
!  ! and finally the velocity vector
!  vel_vec = vel_rad + vel_ang
! END IF
