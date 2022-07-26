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
  rad_vel = model_grid(mod_index)%vel
  vel_vec = rad_vel * pos / norm2(pos)
CASE DEFAULT
 write(*,*) 'velo_vector: the choice of velApprox = ', velApprox, ' is not known...'
 STOP
END SELECT


END SUBROUTINE velo_vector
