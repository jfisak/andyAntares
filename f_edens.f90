! Function to calculate the root of the electron number density
!
! INPUT: model_grid_index(INT): modGrid cell index
!        el_nd(DBL): electron density
! OUTPUT: func(DBL): calculated function
!
SUBROUTINE f_edens(model_grid_index, el_nd, func)


USE types
USE constants

IMPLICIT NONE
INTEGER            :: model_grid_index
INTEGER            :: indexe, indexi, numb_ions
DOUBLE PRECISION   :: temp, el_nd, func, SUMME, SUMMI, frac

DOUBLE PRECISION        :: abundance


temp = model_grid(model_grid_index)%T

SUMME = 0.D0 

DO indexe = 1, n_elements  
 SUMMI = 0.D0
 numb_ions = elements(indexe)%nions
 DO indexi = 1, numb_ions
  CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
  SUMMI = SUMMI + DBLE(indexi - 1) * frac
 END DO
 abundance = model_grid(model_grid_index)%grid_comp(indexe)%abund
 SUMME = SUMME + (abundance / elements(indexe)%atom_mass) * SUMMI
 ! write(*,*) 'f_edens: abundance = ', abundance, ' atmass = ', elements(indexe)%atom_mass
END DO

! we are searching for solution el_nd and because we subtract el_nd 
! write(*,*) 'f_edens: 1 = ', model_grid(model_grid_index)%rho * SUMME, ' el_nd = ', el_nd
func = model_grid(model_grid_index)%rho * SUMME - el_nd  

END SUBROUTINE f_edens
