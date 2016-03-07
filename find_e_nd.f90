SUBROUTINE find_e_nd(model_grid_index, el_nd)
   
! Calculation of electron number density iteratively in 
! model grid cellmodel_grid_index

  USE types

  IMPLICIT NONE
  INTEGER, PARAMETER                 :: max_it = 100
  DOUBLE PRECISION, PARAMETER        :: minc =10.D-05
  INTEGER                            :: loop_index, model_grid_index
  DOUBLE PRECISION                   :: el_nd, el_nd_1, el_nd_2, func1, func2
  DOUBLE PRECISION                   :: diff

  ! Smallest value for the el_nd (can be either 1. or 0.)
  el_nd_1 = 0.D0

!  DO I = 1, n_elements
!     numb_ions = elements(I)%nions
!     e_nd_2 = model_grid(index_model_grid)%rho * (model_grid(model_grid_index)%grid_comp(I)%abund / elements(I)%atom_mass) * (numb_ions - 1)
!  END DO

   ! All H is ionized because and it is larger value for the el_nd 
   el_nd_2 = model_grid(model_grid_index)%rho / mp_g
!   print*, 'hodnota el_nd_2: ', el_nd_2

  ! Debug
  ! print*, el_nd_1, el_nd_2
  ! CALL f_edens(model_grid_index, el_nd_1, func1)
  ! CALL f_edens(model_grid_index, el_nd_2, func2)
  ! print*, func1, func2
  ! STOP

  loop_index = 1
!    print*, "loop_index: ", loop_index
!    print*, "max_it: ", max_it
  DO 
!    print*, "loop_index: ", loop_index
!    print*, "max_it: ", max_it
    IF (loop_index .GT. max_it) THEN
        PRINT*, 'No solution for the electron density found in cell', model_grid_index
        STOP
    END IF
    ! Calculate the function func1 which is the root of the electron number density
    CALL f_edens(model_grid_index, el_nd_1, func1)
    CALL f_edens(model_grid_index, el_nd_2, func2)
!    print*, 'electron densities: el_nd_1 = ', el_nd_1, ', el_nd_2 = ', el_nd_2
    diff = (el_nd_2 - el_nd_1) / (func2 - func1) * func2 
    IF (ABS(diff) .LT. minc) THEN
!       PRINT*, 'Electron number density in cell', model_grid_index, 'equals', el_nd_2
       EXIT
    END IF
    el_nd_1 = el_nd_2
    el_nd_2 = el_nd_2 - diff

    ! Debug
!    print*, 'electron density: ', loop_index, el_nd_1, el_nd_2, '\n', diff
    !print*, 'electron density: ', diff
    !print*, loop_index, func1, func2, diff

    loop_index = loop_index + 1  
   END DO    
  
   el_nd = el_nd_2
  
END SUBROUTINE find_e_nd

