  SUBROUTINE saha_boltzmann_factor(indexe, indexi, temp, sb_factor)
  
  ! Calculate the ratio sb_factor between population of any two exication levels 
  ! indexi and (indexi+1)
  USE types  

  IMPLICIT NONE

  INTEGER                              :: indexe, indexi, model_grid_index
  DOUBLE PRECISION                     :: temp, sb_factor, U1, U2

!  print*, 'Saha fact. called for:',indexe, indexi, temp

  CALL part_fun(indexe, indexi, temp, U1)
  CALL part_fun(indexe, indexi+1, temp, U2)
!  print*, '  Part.func:', U1, U2

  sb_factor = U1 / U2 * saha_const * temp**(-3.D0/2.D0) * &
               EXP( elements(indexe)%ions(indexi)%ion_potential / (BOLK * temp) )
! print*, '  exponenciala: ', EXP( elements(indexe)%ions(indexi)%ion_potential / (BOLK * temp) )
! print*, '  elements(indexe)%ions(indexi%ion_potential = ', elements(indexe)%ions(indexi)%ion_potential, ' BOLK = ', BOLK
! print*, '  temp: ', temp
! print*, '  Ion.pot.', elements(indexe)%ions(indexi)%ion_potential/e_v
! print*, '  Saha Boltz. factor:', sb_factor 

  END SUBROUTINE saha_boltzmann_factor

