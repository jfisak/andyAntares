! Calculate the fraction of atoms of chemical element indexe in ionization
! stage indexi relative to the total number of atoms of that element
! 
! INPUT: indexe(INT): element index
!        indexi(INT): ion index
!        temp(DBL): temperature
!        el_nd(DBL): electron density
! OUTPUT: frac(DBL): calculated ionization fraction
!
SUBROUTINE ionization_fraction(indexe, indexi, temp, el_nd, frac)

USE types 
USE constants

IMPLICIT NONE

INTEGER                                 :: ind_I, ind_J, numb_ions, indexe, indexi
DOUBLE PRECISION                        :: temp, el_nd, frac
DOUBLE PRECISION                        :: tot_N, factor_D, SUMM
DOUBLE PRECISION                        :: sb_factor
LOGICAL                                 :: too_large

! write(*,*)  'START SUBROUTINE ionization_fraction' 
! write(*,*)  'Ion.frac. is called for:',indexe, indexi, temp, el_nd

numb_ions = elements(indexe)%nions
! write(*,*) '  numb.ions:', numb_ions

SELECT CASE(nlte)
! LTE approximation
CASE(0)

 tot_N = 1.D0
 DO ind_I = indexi, numb_ions - 1
  CALL saha_boltzmann_factor(indexe, ind_I, temp, sb_factor, too_large)
  ! write(*,*) 'ionization_fraction: ', ind_I, ' sb_factor = ', sb_factor, tot_N, too_large
  ! IF(too_large .EQV. .TRUE.) CYCLE
  tot_N = tot_N * el_nd * sb_factor
 ! write(*,*) 'ionization_fraction: e  = ', el_nd, ' sf = ', sb_factor, ' tot_N = ', tot_N
 END DO
 
 SUMM = 0.D0
 DO ind_I = 1, numb_ions
  factor_D = 1.D0
  DO ind_J = ind_I, numb_ions - 1
   CALL saha_boltzmann_factor(indexe, ind_J, temp, sb_factor, too_large)
   IF(too_large .EQV. .TRUE.) CYCLE
   ! IF(sb_factor == -1.0) EXIT
   factor_D = factor_D * el_nd * sb_factor
   IF(isnan(factor_D)) THEN
   END IF
  ! write(*,*) 'ionization_fraction: hodnoty: el_nd = ', el_nd, ' ,sb_factor = ', sb_factor, ' , temp = ', temp
  ! write(*,*) 'ionization_fraction: probehl', ind_I, ' a ', ind_J, ' -ty cyklus vypoctu D, factor_D=', factor_D
  END DO
  SUMM = SUMM + factor_D     
  ! write(*,*) 'ionization_fraction: SUMM = ', SUMM
 END DO
 
  !if(SUMM == 0) print*, "ionization_fraction: SUMM = 0..."
   frac = DBLE(tot_N / SUMM)
   ! write(*,*) 'ionization_fraction: fraction = ', frac
   ! IF(frac < 1.D-40) frac = 0
   IF(isnan(frac)) THEN
    write(*,*) 'ionization_fraction: temp = ', temp
    STOP 'ionization_fraction: frac = NaN'
   END IF
! the ionization fractions are read from the input model
CASE(5)
! nothing to be done here
CASE DEFAULT
 write(*,*) 'ionization_fraction: nlte = ', nlte
 STOP 'ionization_fraction: this choice is not possible'
END SELECT

END SUBROUTINE ionization_fraction
