SUBROUTINE freq_from_planck(freq)

  USE types

  IMPLICIT NONE 

  INTEGER                       :: end_loop
  DOUBLE PRECISION, PARAMETER   :: wien_const = 5.879D10 
  DOUBLE PRECISION              :: freq, freq_max, planck, ran_freq, ran_planck, planck_max, ran2
  DOUBLE PRECISION              :: planck_numax, planck_numin
 
! OPEN (UNIT=25, FILE='gauss-3.dat') 
  
  freq_max = wien_const * T_eff

  planck_numin = ( 2.D0 * h * nu_min**3 / light_speed**2  ) * (  1.D0 / ( EXP( (h * nu_min) / (BOLK * T_eff) ) - 1.D0 )  )
  planck_numax = ( 2.D0 * h * nu_max**3 / light_speed**2  ) * (  1.D0 / ( EXP( (h * nu_max) / (BOLK * T_eff) ) - 1.D0 )  )
  planck_max = ( 2.D0 * h * freq_max**3 / light_speed**2  ) * (  1.D0 / ( EXP( (h * freq_max) / (BOLK * T_eff) ) - 1.D0 )  )
 
  IF ((freq_max .GT. nu_min).AND.(freq_max .LT. nu_max)) THEN
     planck_max = planck_max
  ELSE
     
     IF (planck_numax .GT. planck_numin) THEN
        planck_max = planck_numax
     ELSE
        planck_max = planck_numin
     END IF
  END IF
 
  end_loop = 0
  DO WHILE (end_loop .EQ. 0)

     ran_freq = nu_min + (nu_max - nu_min) * ran2(idum)
     ran_planck = ran2(idum) * planck_max
    
     planck =( 2.D0 * h * ran_freq**3 / light_speed**2  ) * (  1.D0 / ( EXP( (h * ran_freq) / (BOLK * T_eff) ) - 1.D0 )  )

     IF ( ran_planck .LT. planck ) THEN
        freq = ran_freq
!        WRITE(25, *) freq, planck
        end_loop = 1
     END IF

  END DO

END SUBROUTINE freq_from_planck
