SUBROUTINE update_grid(iteration)

  ! Calculate electron number density, population number of the ground level
  ! and total population number for every model grid cell for given composition
  ! and corresponding ionization stages.
  USE types

  IMPLICIT NONE    

INTEGER             :: gridcell, indexe, indexi, numb_ions, iteration
DOUBLE PRECISION    :: el_nd, temp
LOGICAL                 :: wasFound

  
write(99,*) 'updating grid'
DO gridcell = 1, n_modelgrid
  IF (model_grid(gridcell)%assoc_cells .GT. 0) THEN
    IF (iteration .EQ. 1) THEN
      ! Calculate electron number density for every model grid cell gridcell
      IF(eldensfile == 0) THEN
       CALL find_e_nd(gridcell, el_nd)
      ELSE
       if(gridcell == 1) CALL read_e_nd()
      END IF
    ELSE
      ! Energy density contribeted to the model grid cell 
       model_grid(gridcell)%J = model_grid(gridcell)%J / model_grid(gridcell)%volume / (4 * pi)
       temp = (model_grid(gridcell)%J * pi / sigma )**(1./4.) 
       ! print*, 'model cell: ', gridcell, ' temperature = ', temp, ' flux = ', model_grid(gridcell)%J, &
       ! ' volume = ', model_grid(gridcell)%volume
       ! print*, 'update_grid: temperature: I = ', I, ' T = ', temp
       model_grid(gridcell)%T = temp
       ! Calculate electron number density for every model grid cell gridcell
       !print*, gridcell, model_grid(gridcell)%J, temp
       CALL find_e_nd(gridcell, el_nd)
       model_grid(gridcell)%J = 0.D0   
    END IF
    IF(eldensfile == 0) THEN
     model_grid(gridcell)%e_dens = el_nd 
     ! write(*,*) 'update_grid: el_nd = ', el_nd
    END IF
    temp = model_grid(gridcell)%T
    !     print*, 'temp and e_nd:', gridcell,  model_grid(gridcell)%rho, temp, el_nd/6.1D14
  ENDIF
END DO
! calculation of population numbers
 IF(iteration == 1) THEN
  CALL find_populations(wasFound)
 END IF
 IF(.NOT. wasFound .OR. iteration > 1) THEN
  DO indexe = 1, n_elements
   numb_ions = elements(indexe)%nions
   DO indexi = 1, numb_ions
    CALL lte_pops(indexe, indexi)
   END DO
  END DO
 END IF
CALL check_pop()
! STOP 'update_grid: testing'
CLOSE(3)

  
END SUBROUTINE update_grid
