 SUBROUTINE angle_aberration(dir1, vel, dir2)

  USE types

  IMPLICIT NONE    

   DOUBLE PRECISION                 :: D_gamma
   DOUBLE PRECISION, DIMENSION(3)   :: dir1, dir2, vel

   D_gamma = 1.D0

   dir2 = dir1 - D_gamma * vel/light_speed * (1.D0 - D_gamma/(D_gamma + 1) * DOT_PRODUCT(dir1,vel)/light_speed)

 END SUBROUTINE angle_aberration
