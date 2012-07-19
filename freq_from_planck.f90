  SUBROUTINE freq_from_planck(freq)

  USE types

  IMPLICIT NONE 

  DOUBLE PRECISION, PARAMETER   :: nu_min= 2.D15, nu_max=3.D15, constant=5.879*10.D10 !nu_min= 1.D14, nu_max=1.D17,
  DOUBLE PRECISION              :: freq, freq_max, planck, ran_freq, ran_planck, planck_max, ran2
  DOUBLE PRECISION              :: planck_numax, planck_numin
    OPEN (UNIT=10, FILE='gauss-3.dat') 
 
10  ran_freq = nu_min + (nu_max - nu_min)*ran2(idum)
    ran_planck = ran2(idum)
    freq_max = constant * T_eff
    planck_numin = ( (2.D0 * h * nu_min**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * nu_min)/(BOLK * T_eff) )-1.D0 )  )
    planck_numax = ( (2.D0 * h * nu_max**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * nu_max)/(BOLK * T_eff) )-1.D0 )  )

     IF ((freq_max .LT. nu_min).AND.(nu_min .LT. freq_max)) THEN
          planck_max = planck_numax
     ELSE
          IF (planck_numax .GT. planck_numin) THEN
              planck_max = planck_numax
          ELSE
              planck_max = planck_numin
          END IF
     END IF

!    planck_max =( (2.D0 * h * freq_max**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * freq_max)/(BOLK * T_eff) )-1.D0 )  )
    
    planck =( (2.D0 * h * ran_freq**3.)/(light_speed**2.)  )*(  1.D0/ ( EXP( (h * ran_freq)/(BOLK * T_eff) )-1.D0 )  )

    IF ( ran_planck .LT. (planck/planck_max) ) THEN
        freq = ran_freq
        WRITE(10, *) freq, planck
    ELSE
        GOTO 10        
    END IF
    

  END SUBROUTINE freq_from_planck
