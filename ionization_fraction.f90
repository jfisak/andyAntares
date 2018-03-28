SUBROUTINE ionization_fraction(indexe, indexi, temp, el_nd, frac)

! Calculate the fraction of atoms of chemical element indexe in ionization
! stage indexi relative to the total number of atoms of that element
USE types 

IMPLICIT NONE

INTEGER                                 :: I, J, numb_ions, indexe, indexi
DOUBLE PRECISION                        :: temp, el_nd, frac
REAL(kind=16)                        :: N, D, SUMM
!DOUBLE PRECISION                        :: N, D, SUMM
DOUBLE PRECISION                        :: sb_factor
LOGICAL                                 :: too_large
DOUBLE PRECISION                        :: divide_factor

!  print*, 'START SUBROUTINE ionization_fraction' 
!  print*, 'Ion.frac. is called for:',indexe, indexi, temp, el_nd

numb_ions = elements(indexe)%nions
!  print*, '  numb.ions:', numb_ions

N = 1.Q0
DO I = indexi, numb_ions - 1
 CALL saha_boltzmann_factor(indexe, I, temp, sb_factor, too_large)
!  print*, 'el_nd: ', el_nd
!   print*, 's-b factor = ', sb_factor
 ! write(*,*) 'saha Boltzman factor: ', I, sb_factor, N, too_large
 IF(too_large .EQV. .TRUE.) CYCLE
 N = N * el_nd * sb_factor
! write(*,*) 'ionization_fraction: e * sf = ', el_nd * sb_factor, ' N = ', N
END DO

SUMM = 0.D0
DO I = 1, numb_ions
 D = 1.Q0
 DO J = I, numb_ions - 1
  CALL saha_boltzmann_factor(indexe, J, temp, sb_factor, too_large)
  IF(too_large .EQV. .TRUE.) CYCLE
  !print*, 'ionization_fraction: sb_factor = ', sb_factor
  ! IF(sb_factor == -1.0) EXIT
  D = D * el_nd * sb_factor
  IF(isnan(D)) THEN
  END IF
! print*, 'vypocet ionization fraction, hodnoty: el_nd = ', el_nd, ' ,sb_factor = ', sb_factor, ' , temp = ', temp
! print*, 'probehl', I, ' a ', J, ' -ty cyklus vypoctu D, D=', D
 END DO
 SUMM = SUMM + D     
! print*, 'hodnota SUMM = ', SUMM
END DO

 !if(SUMM == 0) print*, "ionization_fraction: SUMM = 0..."
  frac = DBLE(N / SUMM)
  IF(frac < 1.D-40) frac = 0
!  print*, '   Ion.frac:', frac
!  print*, 'END SUBROUTINE ionization_fraction' 

END SUBROUTINE ionization_fraction
