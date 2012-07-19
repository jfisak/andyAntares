  SUBROUTINE toy_model() 

  USE types

  IMPLICIT NONE    

  INTEGER                           :: I, J, M
  DOUBLE PRECISION                  :: delta_r, r, delta, delta2

  delta_r = (R_inf - R_star)/n_modelgrid

  DO I=1,n_modelgrid
     r = R_star + I * delta_r
     model_grid(I)%rwind = r
     model_grid(I)%vel   = (r/R_inf)*V_inf
     model_grid(I)%rho   = M_dot / (4.D0 * pi * r**2 * model_grid(I)%vel)     
  END DO


 DO I=1, Ngrid
!    Absolute radius of the grid cell 
     r =SQRT( (cell(I)%corner(1) + cell_width/2.D0)**2 + &
              (cell(I)%corner(2) + cell_width/2.D0)**2 + &
              (cell(I)%corner(3) + cell_width/2.D0)**2)
     delta = 1.D99
     DO J=1,n_modelgrid   
        delta2 = ABS(r -model_grid(J)%rwind)
        IF (delta2 .LT. delta) THEN
            delta = delta2 
            M = J           
        END IF          
     END DO  
     cell(I)%model_index = M     
  END DO
  
  END SUBROUTINE toy_model
