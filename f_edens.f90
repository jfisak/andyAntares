  SUBROUTINE f_edens(model_grid_index, el_nd, func)
  
  ! Function to calculate the root of the electron number density

  USE types

  IMPLICIT NONE
  INTEGER            :: model_grid_index
  INTEGER            :: indexe, indexi, numb_ions
  DOUBLE PRECISION   :: temp, el_nd, func, SUMME, SUMMI, frac

!  print*, 'function FIND_EDENS.F90 start'
!  print*, 'f_edens is called for:', model_grid_index, el_nd

  temp = model_grid(model_grid_index)%T
!  print*, '  cell temp., mass.dens.', temp, model_grid(model_grid_index)%rho

  SUMME = 0.D0 
!  print*, '  numb.elem.', n_elements

  DO indexe = 1, n_elements  
     SUMMI = 0.D0
     numb_ions = elements(indexe)%nions
!     print*, '  numb.ions', indexe, numb_ions
     DO indexi = 1, numb_ions
        CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
        SUMMI = SUMMI + (indexi - 1) * frac
!        print*, '  Ion.frac:', indexi, frac
     END DO
!     print*, model_grid(model_grid_index)%grid_comp(indexe)%abund, elements(indexe)%atom_mass
     SUMME = SUMME + (model_grid(model_grid_index)%grid_comp(indexe)%abund / elements(indexe)%atom_mass) * SUMMI
  END DO

  ! we are searching for solution el_nd and because we subtract el_nd 
!  print*, 'model_grid(model_grid_index)%rho: ', model_grid(model_grid_index)%rho, &
!           'SUMME: ', SUMME, 'el_nd: ', el_nd
  func = model_grid(model_grid_index)%rho * SUMME - el_nd  


!  print*, '  f_edens:', func
!  print*, 'function FIND_EDENS.F90 end'

  END SUBROUTINE f_edens
