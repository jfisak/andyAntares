  SUBROUTINE ionization_fraction(indexe, indexi, temp, el_nd, frac)
  
  ! Calculate the fraction of atoms of chemical element indexe in ionization
  ! stage indexi relative to the total number of atoms of that element
  USE types 

  IMPLICIT NONE

  INTEGER                                 :: I, J, numb_ions, indexe, indexi
  DOUBLE PRECISION                        :: temp, el_nd, frac
  DOUBLE PRECISION                        :: N, D, SUMM, sb_factor

  print*, 'START SUBROUTINE ionization_fraction' 
!  print*, 'Ion.frac. is called for:',indexe, indexi, temp, el_nd

  numb_ions = elements(indexe)%nions
!  print*, '  numb.ions:', numb_ions

  N = 1.D0
  DO I = indexi, numb_ions - 1
     CALL saha_boltzmann_factor(indexe, I, temp, sb_factor)
   print*, 'el_nd: ', el_nd
     N = N * el_nd * sb_factor
     print*, 'saha Boltzman factor: ', I, sb_factor, N
  END DO

  SUMM = 0.D0
  DO I = 1, numb_ions
     D = 1.D0
     DO J = I, numb_ions - 1
        CALL saha_boltzmann_factor(indexe, J, temp, sb_factor)
        D = D * el_nd * sb_factor
!        print*, 'vypocet ionization fraction, hodnoty: el_nd = ', el_nd, ' ,sb_factor = ', sb_factor, ' , temp = ', temp
!        print*, 'probehl', I, ' a ', J, ' -ty cyklus vypoctu D, D=', D
     END DO
     SUMM = SUMM + D     
!      print*, 'hodnota sumy SUMM = ', SUMM
  END DO

  frac = N / SUMM
  print*, '   Ion.frac:', frac
  print*, 'END SUBROUTINE ionization_fraction' 

  END SUBROUTINE ionization_fraction
