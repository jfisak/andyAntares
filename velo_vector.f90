SUBROUTINE velo_vector(pos, mod_index, vel_vec)

USE types
IMPLICIT NONE

DOUBLE PRECISION, DIMENSION(3)                          :: pos
INTEGER                                                 :: mod_index

DOUBLE PRECISION, DIMENSION(3)                          :: vel_vec

DOUBLE PRECISION                                        :: rad_vel


SELECT CASE(model_type)
 ! spherically symmetric models
 CASE(1)
  IF(norm2(pos) < R_star .or. norm2(pos) > R_inf) THEN
   vel_vec = (/ 0.0, 0.0, 0.0 /)
   RETURN
  END IF
  rad_vel = model_grid(mod_index)%vel
  vel_vec = rad_vel * pos / norm2(pos)
 CASE(3)
  IF(mod_index <= n_modelgrid) THEN
   vel_vec = model_grid(mod_index)%vec_vel
  ELSE IF(mod_index == n_modelgrid + 3) THEN ! special treatment for vacuum cells
   vel_vec = (/-2.0, -3.0, -5.0/)
  ELSE
   vel_vec = (/ 0.0, 0.0, 0.0 /)
  END IF
CASE DEFAULT
 write(*,*) 'velo_vector: the choice of velApprox = ', velApprox, ' is not known...'
 STOP
END SELECT

END SUBROUTINE velo_vector
