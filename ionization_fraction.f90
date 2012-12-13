  SUBROUTINE ionization_fraction(indexe, indexi, temp, el_nd, frac)
  
  ! Calculate the fraction of atoms of chemical element indexe in ionization
  ! stage indexi relative to the total number of atoms of that element
  USE types 

  IMPLICIT NONE

  INTEGER                                 :: I, J, numb_ions, indexe, indexi
  DOUBLE PRECISION                        :: temp, el_nd, frac
  DOUBLE PRECISION                        :: N, D, SUMM, sb_factor

!  print*, 'Ion.frac. is called for:',indexe, indexi, temp, el_nd

  numb_ions = elements(indexe)%nions
!  print*, '  numb.ions:', numb_ions

  N = 1.D0
  DO I = indexi, numb_ions - 1
     CALL saha_boltzmann_factor(indexe, I, temp, sb_factor)
     N = N * el_nd * sb_factor
!     print*, '   ', I, sb_factor, N
  END DO

  SUMM = 0.D0
  DO I = 1, numb_ions
     D = 1.D0
     DO J = I, numb_ions - 1
        CALL saha_boltzmann_factor(indexe, J, temp, sb_factor)
        D = D * el_nd * sb_factor
     END DO
     SUMM = SUMM + D     
  END DO

  frac = N / SUMM
!  print*, '   Ion.frac:', frac
   

  END SUBROUTINE ionization_fraction
