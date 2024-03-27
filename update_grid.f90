SUBROUTINE update_grid(iteration)

! Calculate electron number density, population number of the ground level
! and total population number for every model grid cell for given composition
! and corresponding ionization stages.
USE types
USE constants

  IMPLICIT NONE    

INTEGER             :: indexe, indexi, numb_ions, iteration
DOUBLE PRECISION    :: el_nd, temp
LOGICAL                 :: wasFound

INTEGER                                 :: N_single, N_zbytek
INTEGER                                 :: cur_mgi
INTEGER                                 :: my_start, my_end

LOGICAL                                 :: propmod_file_exists
CHARACTER(60)                           :: propmod_file

INTEGER                                 :: cur_n_assoccells

DOUBLE PRECISION, DIMENSION(n_modelgrid + add_mg)    :: cur_j, cur_temp, cur_elnd

cur_j(:) = 0.D0
cur_temp(:) = 0.D0
cur_elnd(:) = 0.D0
  
write(99,*) 'updating grid'
#if mpi == 1
 N_single = n_modelgrid/n_tasks
 N_zbytek = n_modelgrid - n_tasks * N_single
 IF(my_rank <= N_zbytek - 1) THEN
  my_start = my_rank * (N_single + 1) + 1
  my_end = my_rank * (N_single + 1) + N_single
 ELSE IF(N_zbytek == 0) THEN
  my_start = my_rank * (N_single + 1) + 1
  my_end = my_rank * (N_single + 1) + N_single
 ELSE
  my_start = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + 1
  my_end = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + N_single +1
 END IF
#else
 my_start = 1
 my_end = n_modelgrid
#endif

write(*,*) 'update_grid: my_start = ', my_start, ' my_end = ', my_end

DO cur_mgi = my_start, my_end
 cur_n_assoccells = model_grid(cur_mgi)%assoc_cells
 IF (cur_n_assoccells .GT. 0) THEN
  IF (iteration == 1) THEN
   ! Calculate electron number density for every model grid cell gridcell
   IF(eldensfile == 0) THEN
    CALL find_e_nd(cur_mgi, el_nd)
   END IF
   ! reading the temperature structure
   cur_temp(cur_mgi) = model_grid(cur_mgi)%T
   ! write(*,*) 'update_grid: cur_temp = ', cur_temp(cur_mgi)
  ELSE ! iteration > 1
   ! Energy density contribeted to the model grid cell 
   cur_j(cur_mgi) = model_grid(cur_mgi)%J / model_grid(cur_mgi)%volume / (4 * pi)
   temp = (model_grid(cur_mgi)%J * pi / sigma )**(1./4.) 
   cur_temp(cur_mgi)  = temp
   ! Calculate electron number density for every model grid cell gridcell
   CALL find_e_nd(cur_mgi, el_nd)
   model_grid(cur_mgi)%J = 0.D0   
  END IF ! test for the first iteration
  IF(eldensfile == 0) THEN
   cur_elnd(cur_mgi) = el_nd 
  END IF
  write(propmod_file,"(A, A12)") TRIM(outputfolder), '/propmod.dat'
  INQUIRE(FILE=propmod_file, EXIST=propmod_file_exists)

  ! write(*,*) 'update_grid: saved_grid = ', saved_grid, ' propmod_file_exists = ', propmod_file_exists
  IF(saved_grid == 1 .and. propmod_file_exists) THEN
  ! nothing
  ELSE
  ! if we calculate the condition only from electron density
   ! write(*,*) 'update_grid: enable_diffusion = ', enable_diffusion
   IF(enable_diffusion == 1) THEN
    CALL diffusion_approximation(cur_mgi)
   END IF
  END IF
 END IF ! cur_n_assoccells > 0
END DO

#if mpi == 1
 CALL mpi_reduce(cur_j(:), model_grid(:)%j, n_modelgrid + add_mg, &
  & mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
 CALL MPI_REDUCE(cur_temp(:), model_grid(:)%T, n_modelgrid + add_mg, &
  & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_REDUCE(cur_elnd(:), model_grid(:)%e_dens, n_modelgrid + add_mg, &
  & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 ! distribute all the data among all of the processes
 CALL MPI_BCAST(model_grid(:)%e_dens, n_modelgrid + add_mg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_BCAST(model_grid(:)%j, n_modelgrid + add_mg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_BCAST(model_grid(:)%T, n_modelgrid + add_mg, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
#endif


! if(my_rank == 0) then
!  write(*,*) 'update_grid: cur_e_dens = ', model_grid(:)%e_dens
! end if

! STOP 'update_grid: testing'
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
