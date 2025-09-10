! this sbr calculates a velocity vector saved in the current model grid
!
! input: pos -- DBLE(3) vector of a position
!        mod_index -- INT modGrid cell index
! output: vel_vec -- DBLE(3) vector of a velocity
!
! RETURN point 2X
SUBROUTINE velo_vector(pos, mod_index, vel_vec)

USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: pos
INTEGER                                                 :: mod_index

INTEGER                                                 :: cur_pgi

DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: vel_vec

DOUBLE PRECISION                                        :: rad_vel

DOUBLE PRECISION                                        :: vr, vtheta
DOUBLE PRECISION                                        :: coor_x, coor_y, coor_z
DOUBLE PRECISION                                        :: cosphi, sinphi, costheta, sintheta
DOUBLE PRECISION                                        :: radius, phi, theta
DOUBLE PRECISION, DIMENSION(const_dimofspace)           :: cur_n

! calculation of the model index if it is not defined in the input
IF(mod_index == no_modcell) THEN
 CALL find_dyn_cell1(pos, cur_pgi)
 mod_index = dyn_cell(cur_pgi)%model_index
END IF

SELECT CASE(model_type)
 ! spherically symmetric models
 CASE(1)
  IF(norm2(pos) <= R_star) THEN
   vel_vec = (/ 0.D0, 0.D0, 0.D0 /)
   IF(velApprox == 3) THEN
    cur_n = pos/norm2(pos)
    vel_vec = V_star * cur_n
    RETURN
   END IF
   IF(inputmodel == 5) THEN
    vel_vec = V_star * pos / norm2(pos)
   END IF
   RETURN
  ELSE IF(norm2(pos) > R_inf) THEN
   vel_vec = V_inf * pos / norm2(pos)
   IF(velApprox == 1) THEN
    rad_vel = V_inf * (1.0 - R_star / R_inf)**beta
    vel_vec = rad_vel * pos/NORM2(pos)
   END IF
   RETURN
  END IF
  IF(inputmodel /= 5) THEN
   rad_vel = model_grid(mod_index)%vel
   vel_vec = rad_vel * pos / norm2(pos)
  ELSE IF(inputmodel == 5) THEN
   vr = model_grid(mod_index)%vel
   vtheta = model_grid(mod_index)%velang

   coor_x = pos(ind_x)
   coor_y = pos(ind_y)
   coor_z = pos(ind_z)

   radius = norm2(pos)
   theta = acos(coor_z/radius)
   phi = atan2(coor_y,coor_x)

   sintheta = sin(theta)
   costheta = cos(theta)
   sinphi = sin(phi)
   cosphi = cos(phi)

   vel_vec(ind_x) = vr * sintheta * cosphi + vtheta * costheta * cosphi 
   vel_vec(ind_y) = vr * sintheta * sinphi + vtheta * costheta * sinphi
   vel_vec(ind_z) = vr * costheta - vtheta * sintheta
  END IF
 !________________________________________________________
 ! 2D model with radial velocity and tangential velocity
 CASE(2)
  vr = model_grid(mod_index)%vel
  vtheta = model_grid(mod_index)%velang

  coor_x = pos(ind_x)
  coor_y = pos(ind_y)
  coor_z = pos(ind_z)

  ! cosphi = coor_y / norm2(pos)
  ! sinphi = coor_x / norm2(pos)
  
  radius = norm2(pos)
  theta = acos(coor_z/radius)
  phi = atan2(coor_y,coor_x)

  sintheta = sin(theta)
  costheta = cos(theta)
  sinphi = sin(phi)
  cosphi = cos(phi)

  ! costheta = coor_z/norm2(pos)
  ! sintheta = sqrt(coor_x**2+coor_y**2)/norm2(pos)

  ! output vector
  ! vel_vec(1) = vr * costheta * cosphi - vtheta * sintheta * cosphi
  ! vel_vec(2) = vr * costheta * sinphi - vtheta * sintheta * sinphi
  ! vel_vec(3) = vr * sintheta + vtheta * costheta

  ! vel_vec(ind_x) = vr * costheta * cosphi - vtheta * cosphi
  vel_vec(ind_x) = vr * sintheta * cosphi + vtheta * costheta * cosphi 
  vel_vec(ind_y) = vr * sintheta * sinphi + vtheta * costheta * sinphi
  vel_vec(ind_z) = vr * costheta - vtheta * sintheta


  ! vel_vec(ind_y) = vr * costheta * sinphi - vtheta * sinphi
  ! vel_vec(ind_z) = vr * sintheta
 
 !________________________________________________________
 CASE(3)
  IF(mod_index <= n_modelgrid) THEN
   vel_vec = model_grid(mod_index)%vec_vel
  ELSE IF(mod_index == photosphere_index) THEN
   vel_vec = (/ 0.0, 0.0, 0.0 /)
  ELSE IF(mod_index == outerspace_index) THEN
   vel_vec = R_inf * pos/norm2(pos)
  ELSE IF(mod_index == vacuum_index) THEN ! special treatment for vacuum cells
   vel_vec = (/-2.0, -3.0, -5.0/)
  END IF
 !________________________________________________________
 ! 1D model with radial velocity and tangential velocity
 CASE(6)
  vr = model_grid(mod_index)%vel
  vtheta = model_grid(mod_index)%velang

  coor_x = pos(ind_x)
  coor_y = pos(ind_y)
  coor_z = pos(ind_z)

  ! cosphi = coor_y / norm2(pos)
  ! sinphi = coor_x / norm2(pos)
  
  radius = norm2(pos)
  theta = acos(coor_z/radius)
  phi = atan2(coor_y,coor_x)

  sintheta = sin(theta)
  costheta = cos(theta)
  sinphi = sin(phi)
  cosphi = cos(phi)


  ! vel_vec(ind_x) = vr * costheta * cosphi - vtheta * cosphi
  vel_vec(ind_x) = vr * sintheta * cosphi + vtheta * costheta * cosphi 
  vel_vec(ind_y) = vr * sintheta * sinphi + vtheta * costheta * sinphi
  vel_vec(ind_z) = vr * costheta - vtheta * sintheta


  ! vel_vec(ind_y) = vr * costheta * sinphi - vtheta * sinphi
  ! vel_vec(ind_z) = vr * sintheta
 
 !________________________________________________________
 CASE DEFAULT
 write(*,*) 'velo_vector: the choice of model_type = ', model_type, ' is not known...'
 STOP
END SELECT

END SUBROUTINE velo_vector
