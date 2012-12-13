SUBROUTINE do_spectrum(n_pack)

  USE types
  
  IMPLICIT NONE    

  INTEGER                               :: I, n_pack, pack_index, nubin
  DOUBLE PRECISION                      :: delta_nu, delta_e, freq, lambda, flambda, planck, ls_A, frequency
  TYPE(spec_type), DIMENSION(n_nubin)   :: spectrum

  OPEN (UNIT=19, FILE='spec.dat')     
 
  ! Set up the frequency grid to extract spectrum
  print*, 'SETUP FREQ GRID'
  delta_nu = (nu_max - nu_min / 2.D0) / n_nubin
  DO I= 1, n_nubin 
     spectrum(I)%freq = nu_min + (I - 1) * delta_nu
     spectrum(I)%flux = 0.D0           
     spectrum(I)%esc = 0
     !print*, i, spectrum(i)%freq, spectrum(i)%flux
  END DO
    
  ! Loop over all packets
  print*, 'BIN PACKETS' 
  DO pack_index = 1, n_pack
     ! And take all which actually escaped
     IF (package(pack_index)%typ .EQ. type_escaped) THEN
        freq = package(pack_index)%freq_rf
        ! Only bin those packets which are in the allowed frequency range
        IF ((freq .GT. nu_min) .AND. (freq .LT. nu_max)) THEN
           nubin = floor( (freq - nu_min) / delta_nu ) + 1
           delta_e = (package(pack_index)%e_rf / delta_nu) / (4.D0 * pi * (100.D0 * parsec)**2) !put the star to 100 parsecs
           spectrum(nubin)%flux = spectrum(nubin)%flux + delta_e
           spectrum(nubin)%esc = spectrum(nubin)%esc + 1
        ENDIF
     END IF
  END DO


! DO I = 1, n_nubin
!     frequency = spectrum(I)%freq
!     planck =( 2.D0 * h * frequency**3 / light_speed**2  ) * (  1.D0 / ( EXP( (h * frequency) / (BOLK * T_eff) ) - 1.D0 )  )
!     WRITE(19,*)  frequency, spectrum(I)%flux, spectrum(I)%esc, planck
! END DO

  print*, 'WRITE TO FILE'
  
  ls_A = light_speed * 1.D8

  DO I = 1, n_nubin  
     lambda = ls_A / spectrum(I)%freq
     flambda = spectrum(I)%flux * ( ls_A / lambda**2 )
     planck =( 2.D0 * h * ls_A**2 / lambda**5 ) * &
             ( 1.D0 / ( EXP( h * ls_A / (lambda * BOLK * T_eff) ) - 1.D0) )
     WRITE(19,*) lambda, flambda, planck, flambda/planck, spectrum(I)%esc
    
  END DO


END SUBROUTINE do_spectrum
