SUBROUTINE update_grid(iteration)

  ! Calculate electron number density, population number of the ground level
  ! and total population number for every model grid cell for given composition
  ! and corresponding ionization stages.
  USE types

  IMPLICIT NONE    

INTEGER             :: gridcell, indexe, indexi, numb_ions, iteration
DOUBLE PRECISION    :: el_nd, temp
LOGICAL                 :: wasFound

LOGICAL                                 :: propmod_file_exists
CHARACTER(60)                           :: propmod_file
  
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
     model_grid(gridcell)%T = temp
     ! Calculate electron number density for every model grid cell gridcell
     CALL find_e_nd(gridcell, el_nd)
     model_grid(gridcell)%J = 0.D0   
  END IF
  IF(eldensfile == 0) THEN
   model_grid(gridcell)%e_dens = el_nd 
   ! write(*,*) 'update_grid: el_nd = ', el_nd
  END IF
  write(propmod_file,"(A, A12)") TRIM(outputfolder), '/propmod.dat'
  INQUIRE(FILE=propmod_file, EXIST=propmod_file_exists)

  IF(saved_grid == 1 .and. propmod_file_exists) THEN
  ! nothing
  ELSE
  ! if we calculate the condition only from electron density
   IF(enable_diffusion == 1) THEN
   CALL diffusion_approximation(gridcell)
   END IF
  END IF
  temp = model_grid(gridcell)%T
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
! CALL save_rates()
! STOP 'update_grid: testing'
CLOSE(3)

  
END SUBROUTINE update_grid
