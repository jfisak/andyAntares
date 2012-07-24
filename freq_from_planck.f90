SUBROUTINE freq_from_planck(freq)

  USE types

  IMPLICIT NONE 

  INTEGER                       :: end_loop
  DOUBLE PRECISION, PARAMETER   :: constant=5.879*10.D10 
  DOUBLE PRECISION              :: freq, freq_max, planck, ran_freq, ran_planck, planck_max, ran2
  DOUBLE PRECISION              :: planck_numax, planck_numin
  OPEN (UNIT=10, FILE='gauss-3.dat') 

  freq_max = constant * T_eff
  planck_numin = ( (2.D0 * h * nu_min**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * nu_min)/(BOLK * T_eff) )-1.D0 )  )
  planck_numax = ( (2.D0 * h * nu_max**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * nu_max)/(BOLK * T_eff) )-1.D0 )  )

  IF ((freq_max .GT. nu_min).AND.(freq_max .LT. nu_max)) THEN
     planck_max = planck_numax
  ELSE
     IF (planck_numax .GT. planck_numin) THEN
        planck_max = planck_numax
     ELSE
        planck_max = planck_numin
     END IF
  END IF
 
  end_loop=0
  DO WHILE (end_loop .EQ. 0)
     ran_freq = nu_min + (nu_max - nu_min)*ran2(idum)
     ran_planck = ran2(idum)
    
     planck =( (2.D0 * h * ran_freq**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * ran_freq)/(BOLK * T_eff) )-1.D0 )  )

     IF ( ran_planck .LT. (planck/planck_max) ) THEN
        freq = ran_freq
        end_loop = 1
        WRITE(10, *) freq, planck
     ELSE
     END IF
  END DO

END SUBROUTINE freq_from_planck
