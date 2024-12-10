SUBROUTINE find_e_nd(model_grid_index, el_nd)
   
! Calculation of electron number density iteratively in 
! model grid cellmodel_grid_index

  USE types
USE constants

  IMPLICIT NONE
  INTEGER, PARAMETER                 :: max_it = 100
  DOUBLE PRECISION, PARAMETER        :: minc =10.D-05
  INTEGER                            :: loop_index, model_grid_index
  DOUBLE PRECISION                   :: el_nd, el_nd_1, el_nd_2, func1, func2
  DOUBLE PRECISION                   :: diff

  IF(model_grid(model_grid_index)%rho == 0.00) THEN
   el_nd = 0.00
   write(*,*) 'find_e_nd: el_nd set to zero'
   RETURN
  END IF

  ! Smallest value for the el_nd (can be either 1. or 0.)
  el_nd_1 = 1.D0


   ! All H is ionized because and it is larger value for the el_nd 
   el_nd_2 = model_grid(model_grid_index)%rho / const_mp_g
   ! write(*,*) 'find_e_nd: el_nd_2 = ', el_nd_2

  

  loop_index = 1
  DO 
    IF (loop_index .GT. max_it) THEN
        write(*,*) 'find_e_nd: No solution for the electron density found in cell', model_grid_index
        write(*,*) 'temperature = ', model_grid(model_grid_index)%T
        STOP
    END IF
    ! Calculate the function func1 which is the root of the electron number density
    ! write(*,*) 'find_e_nd: calling f_edens for el_nd_1 = ', el_nd_1, ' el_nd_2 = ', el_nd_2
    CALL f_edens(model_grid_index, el_nd_1, func1)
    CALL f_edens(model_grid_index, el_nd_2, func2)
    ! write(*,*) , 'electron densities: func1 = ', func1, ', func2 = ', func2
    diff = (el_nd_2 - el_nd_1) / (func2 - func1) * func2 
    ! write(*,*) 'find_e_nd: el_nd_1 = ', el_nd_1, ' el_nd_2 = ', el_nd_2
    IF (ABS(el_nd_1 / el_nd_2 - 1.D0) .LT. minc) THEN
       EXIT
    END IF
    el_nd_1 = el_nd_2
    el_nd_2 = el_nd_2 - diff

    ! Debug
    ! write(*,*)  'electron density: ', loop_index, el_nd_1, el_nd_2
    ! write(*,*)  'electron density: diff = ', diff
    ! write(*,*)  loop_index, func1, func2, diff

    loop_index = loop_index + 1  
   END DO    
  
   el_nd = el_nd_2
  
END SUBROUTINE find_e_nd

