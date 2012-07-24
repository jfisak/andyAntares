SUBROUTINE angle_aberration(dir1, vel, dir2)

  ! Take a rest-frame direction dir1 and velocity vel and calculate
  ! the angle-aberration to the cmf direction dir2, as described
  ! by Mihalas & Mihalas Eq. 89.6
  ! CMF to RF trafo requires to pass a negative velocity Eq. 89.8
 
  USE types
  
  IMPLICIT NONE    

  DOUBLE PRECISION                 :: D_gamma
  DOUBLE PRECISION, DIMENSION(3)   :: dir1, dir2, vel
  
  D_gamma = 1.D0
   
  ! Formula as given by Mihalas & Mihalas, accurate to (v/c)**2
  dir2 = (dir1 - D_gamma * vel/light_speed * &
       (1.D0 - D_gamma/(D_gamma + 1) * DOT_PRODUCT(dir1,vel)/light_speed)) / &
       (D_gamma*(1-DOT_PRODUCT(dir1,vel)/light_speed))


  ! Simplified version accurate to v/c as given by Abott & Luy and
  ! Mazzali & Lucy. Use this for the comparison to my 1D code.
  dir2 = (dir1 - D_gamma * vel/light_speed) / (D_gamma*(1-DOT_PRODUCT(dir1,vel)/light_speed))



 END SUBROUTINE angle_aberration
