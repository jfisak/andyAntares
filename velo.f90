
SUBROUTINE velo(pack_index,vel_vec)

 USE types

 IMPLICIT NONE    

 INTEGER                           :: pack_index
 DOUBLE PRECISION                  :: vel_radial, vec_length
 DOUBLE PRECISION, DIMENSION(3)    :: vel_vec
 ! Petr Kurfurst's disk model variables
 INTEGER                           :: pack_mi, get_package_model_index
 DOUBLE PRECISION, DIMENSION(3)    :: vel_rad, vel_ang, phot_pos
 DOUBLE PRECISION                  :: vel_rad_norm, vel_ang_norm

 SELECT CASE(velApprox)
 ! homologous expansion
 CASE(0)
  ! vel_radial = V_inf/R_inf * vec_length(package(pack_index)%pos)
  IF(vec_length(package(pack_index)%pos) > R_inf) THEN
   vel_vec = (/ 0.0, 0.0, 0.0/)
   RETURN
  END IF
  vel_radial = V_inf/R_inf * vec_length(package(pack_index)%pos)
  vel_vec = package(pack_index)%pos/vec_length(package(pack_index)%pos) * vel_radial
 ! the beta velocity law
 CASE(1)
  vel_radial = V_inf * (1.D0 - R_star / norm2(package(pack_index)%pos))**beta
  vel_vec = package(pack_index)%pos/vec_length(package(pack_index)%pos) * vel_radial
  ! write(*,*) 'velo: pos = ', package(pack_index)%pos 
 CASE DEFAULT
  STOP 'velo: this velocity structure is not known'
 END SELECT
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! petr kurfurst's disk model
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF ((model_type .EQ. 2) .AND. (inputmodel .EQ. 1)) THEN
  ! we have to know the velocity of matter in the given point
  pack_mi = get_package_model_index(pack_index)
  vel_rad_norm = model_grid(pack_mi)%vel
  vel_ang_norm = model_grid(pack_mi)%velang
  ! the photon position
  phot_pos = package(pack_index)%pos
  ! now we compute given vectors
  vel_rad = (/ vel_rad_norm * phot_pos(1) / vec_length(phot_pos(1)), &
             vel_rad_norm * phot_pos(2) / vec_length(phot_pos(1)), 0.D0 /)
  vel_ang = (/ - vel_ang_norm * phot_pos(2) / vec_length(phot_pos(1)), &
             vel_ang_norm * phot_pos(1) / vec_length(phot_pos(1)), 0.D0 /)
  ! and finally the velocity vector
  vel_vec = vel_rad + vel_ang
 END IF
 IF(norm2(vel_vec) > light_speed) THEN
  write(*,*) 'velo: pack_index = ', pack_index, ' position = ', vec_length(package(pack_index)%pos)/R_sun
  STOP 'velocity is larger than the speed of light'
 END IF
END SUBROUTINE velo
