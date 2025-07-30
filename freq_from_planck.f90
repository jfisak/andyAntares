! generates a random frequency based on plank distribution for a temperature T
!
! INPUT: temperature(DBL): radiation temperature
! OUTPUT: freq(DBL): randomly generated frequency
!
SUBROUTINE freq_from_planck(freq, temperature)

USE types
USE constants

IMPLICIT NONE 

DOUBLE PRECISION                      :: temperature
INTEGER                               :: end_loop
DOUBLE PRECISION, PARAMETER           :: wien_const = 5.879D10 
DOUBLE PRECISION                      :: freq, freq_max, planck, ran_freq, ran_planck, planck_max, ran2
DOUBLE PRECISION                      :: planck_numax, planck_numin
DOUBLE PRECISION                      :: wale_start, wale_end 
DOUBLE PRECISION                      :: nu_max, nu_min
DOUBLE PRECISION                      :: lfreq
DOUBLE PRECISION, PARAMETER           :: delta_f = 20 ! Angstrom


IF(oneline) THEN
 lfreq = linelist(1)%freq
 wale_start = lfreq - delta_f
 wale_end = lfreq + delta_f
ELSE
 wale_start = 200   ! in Angstroms
 wale_end = 20000   ! in Angstroms
END IF

nu_max = const_c / (wale_start * 1.D-8)
nu_min = const_c / (wale_end * 1.D-8)
! write(*,*) 'freq_from_planck: nu_max = ', nu_max, ' nu_min = ', nu_min
 
! OPEN (UNIT=25, FILE='gauss-3.dat') 
  
freq_max = wien_const * temperature

planck_numin = ( 2.D0 * const_h * nu_min**3 / const_c**2  ) &
 * (  1.D0 / ( EXP( (const_h * nu_min) / (const_kB * temperature) ) - 1.D0 ))
planck_numax = ( 2.D0 * const_h * nu_max**3 / const_c**2  ) &
 * (  1.D0 / ( EXP( (const_h * nu_max) / (const_kB * temperature) ) - 1.D0 ))
planck_max = ( 2.D0 * const_h * freq_max**3 / const_c**2  ) &
 * (  1.D0 / ( EXP( (const_h * freq_max) / (const_kB * temperature) ) - 1.D0 ))

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

 planck =( 2.D0 * const_h * ran_freq**3 / const_c**2  ) &
  * (  1.D0 / ( EXP( (const_h * ran_freq) / (const_kB * temperature) ) - 1.D0 )  )

 IF ( ran_planck .LT. planck ) THEN
    freq = ran_freq
     WRITE(*, *) freq, planck
    end_loop = 1
 END IF

END DO
!CLOSE(25)
END SUBROUTINE freq_from_planck
