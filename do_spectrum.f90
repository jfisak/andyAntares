 SUBROUTINE do_spectrum(n_pack)

  USE types

  IMPLICIT NONE    

    INTEGER                               :: i, n_pack, pack_index, nubin
    DOUBLE PRECISION                      :: delta_nu, delta_e, lambda, planck
    TYPE(spec_type), DIMENSION(n_nubin)   :: spectrum

    OPEN (UNIT=9, FILE='spec.dat')     
 
    !Set up the frequency grid to extract spectrum
    print*, 'SETUP FREQ GRID'
    delta_nu = (nu_max - nu_min) / n_nubin
    DO i=1,n_nubin 
      spectrum(i)%freq = nu_min+(i-1)*delta_nu
      spectrum(i)%flux = 0.D0           
    END DO
    
    !Loop over all packets
    print*, 'BIN PACKETS' 
    DO pack_index=1,n_pack
       !And take all which actually escaped
       IF (package(pack_index)%typ .EQ. type_escaped) THEN
           nubin=floor((package(pack_index)%freq_rf-nu_min)/delta_nu)+1
	   delta_e= (package(pack_index)%e_rf / delta_nu) / (4*pi*(100*parsec)**2) !put the star to 100 parsecs
           spectrum(nubin)%flux = spectrum(nubin)%flux + delta_e  
       END IF
    END DO


    print*, 'WRITE TO FILE'
    DO i=1,n_nubin  
      lambda = (light_speed/spectrum(i)%freq)*1.D8
      planck =( (2.D0 * h * (light_speed*1.D8)**2)/(lambda**5)  )* &
              (1.D0/ ( EXP( (h * light_speed * 1.D8)/(lambda * BOLK * T_eff) )-1.D0 )  ) 
      WRITE(9,*) lambda, spectrum(i)%flux*((light_speed*1.D8)/lambda**2), planck, &
                 spectrum(i)%flux*((light_speed*1.D8)/lambda**2)/planck 
!      WRITE(9,*) spectrum(i)%freq, spectrum(i)%flux      
    END DO



 END  SUBROUTINE do_spectrum
