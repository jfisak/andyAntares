SUBROUTINE velo_vector(pos, mod_index, vel_vec)

USE types
USE constants
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: pos
INTEGER                                                 :: mod_index

DOUBLE PRECISION, DIMENSION(3)                          :: vel_vec

DOUBLE PRECISION                                        :: rad_vel

DOUBLE PRECISION                                        :: vr, vtheta
DOUBLE PRECISION                                        :: x, y, z
DOUBLE PRECISION                                        :: cosphi, sinphi, costheta, sintheta


SELECT CASE(model_type)
 ! spherically symmetric models
 CASE(1)
  IF(norm2(pos) < R_star .or. norm2(pos) > R_inf) THEN
   vel_vec = (/ 0.0, 0.0, 0.0 /)
   RETURN
  END IF
  rad_vel = model_grid(mod_index)%vel
  vel_vec = rad_vel * pos / norm2(pos)
  ! write(37,*) norm2(pos)/R_star, rad_vel, mod_index
 !________________________________________________________
 ! 2D model with radial velocity and tangential velocity
 CASE(2)
  vr = model_grid(mod_index)%vel
  vtheta = model_grid(mod_index)%velang

  x = pos(1)
  y = pos(2)
  z = pos(3)

  cosphi = y / norm2(pos)
  sinphi = x / norm2(pos)
  
  costheta = sqrt(x**2+y**2)/norm2(pos)
  sintheta = z/norm2(pos)

  ! output vector
  vel_vec(1) = vr * costheta * cosphi - vtheta * sintheta * cosphi
  vel_vec(2) = vr * costheta * sinphi - vtheta * sintheta * sinphi
  vel_vec(3) = vr * sintheta + vtheta * costheta
 
 !________________________________________________________
 CASE(3)
  IF(mod_index <= n_modelgrid) THEN
   vel_vec = model_grid(mod_index)%vec_vel
  ELSE IF(mod_index == n_modelgrid + 3) THEN ! special treatment for vacuum cells
   vel_vec = (/-2.0, -3.0, -5.0/)
  ELSE
   vel_vec = (/ 0.0, 0.0, 0.0 /)
  END IF
 CASE DEFAULT
 write(*,*) 'velo_vector: the choice of model_type = ', model_type, ' is not known...'
 STOP
END SELECT

END SUBROUTINE velo_vector
