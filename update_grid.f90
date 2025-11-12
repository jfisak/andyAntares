SUBROUTINE update_grid(iteration)

! Calculate electron number density, population number of the ground level
! and total population number for every model grid cell for given composition
! and corresponding ionization stages.
USE MPI
USE types
USE constants

  IMPLICIT NONE    

INTEGER             :: indexe, indexi, numb_ions, iteration
DOUBLE PRECISION    :: el_nd, temp
LOGICAL                 :: wasFound

INTEGER                                 :: N_single, N_zbytek, N_tot_zbytek
INTEGER                                 :: cur_mgi
INTEGER                                 :: my_start, my_end
INTEGER                                 :: status(MPI_STATUS_SIZE)

LOGICAL                                 :: propmod_file_exists
CHARACTER(filename_length)                           :: propmod_file
INTEGER                                 :: ind_I

INTEGER                                 :: cur_n_assoccells
! DOUBLE PRECISION                        :: test_temp

DOUBLE PRECISION, DIMENSION(n_modelgrid + add_mg)    :: cur_j, cur_temp, cur_e_dens
DOUBLE PRECISION, DIMENSION(n_modelgrid + add_mg)    :: recv_j, recv_temp, recv_e_dens

cur_j(:) = 0.D0
recv_j(:) = 0.D0
cur_temp(:) = 0.D0
recv_temp(:) = 0.D0
cur_e_dens(:) = 0.D0
recv_e_dens(:) = 0.D0
  
write(99,*) 'updating grid'
#if mpi == 1
 N_single = (n_modelgrid)/n_tasks
 N_zbytek = n_modelgrid - n_tasks * N_single
 IF(N_zbytek /= 0) N_tot_zbytek = (N_zbytek + 1) * (N_single + 1) + N_zbytek
 write(*,*) 'update_grid: N_single = ', N_single, ' N_zbytek = ', N_zbytek
 IF(my_rank <= N_zbytek - 1) THEN
  my_start = my_rank * N_single  + my_rank + 1
  my_end = (my_rank + 1) * N_single + my_rank + 1
 ELSE IF(N_zbytek == 0) THEN
  my_start = my_rank * N_single + 1
  my_end = (my_rank + 1) * N_single
 ELSE IF(my_rank > N_zbytek - 1) THEN
  my_start = my_rank * N_single  + N_zbytek + 1
  my_end = (my_rank + 1) * N_single + N_zbytek
 END IF
 IF(my_rank == n_tasks - 1) THEN
  my_end = n_modelgrid
 END IF
#else
 my_start = 1
 my_end = n_modelgrid
#endif

write(*,*) 'update_grid: my_rank = ', my_rank, ' n_modelgrid = ', n_modelgrid
write(*,*) 'update_grid: my_start = ', my_start, ' my_end = ', my_end
CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

DO cur_mgi = my_start, my_end
 cur_n_assoccells = model_grid(cur_mgi)%assoc_cells
 ! write(*,*) 'update_grid: my_rank = ', my_rank, ' cur_n_assoccells = ', cur_n_assoccells
 IF (cur_n_assoccells > 0 .and. cur_mgi <= n_modelgrid) THEN
  ! write(*,*) 'update_grid: temp = ', model_grid(cur_mgi)%T
  IF (iteration == 1) THEN
   ! Calculate electron number density for every model grid cell gridcell
   IF(eldensfile == 0) THEN
    ! write(*,*) 'update_grid: cur_mgi = ', cur_mgi, ' n_modelgrid = ', n_modelgrid
    ! write(*,*) 'update_grid: temp = ', model_grid(cur_mgi)%T
    IF(nlte /= 5) THEN ! nlte = 5 == the electron densities are read from the model file
     CALL find_e_nd(cur_mgi, el_nd)
     cur_e_dens(cur_mgi) = el_nd
    END IF
   END IF
   ! reading the temperature structure
   cur_temp(cur_mgi) = model_grid(cur_mgi)%T
  !___________________________________________________________________________________________
  ! iteration > 1
  ELSE 
   ! Energy density contributed to the model grid cell 
   cur_j(cur_mgi) = model_grid(cur_mgi)%J / model_grid(cur_mgi)%volume / (4 * const_pi)
   temp = (model_grid(cur_mgi)%J * const_pi / const_stefbolz )**(1./4.) 
   cur_temp(cur_mgi)  = temp
   ! Calculate electron number density for every model grid cell gridcell
   CALL find_e_nd(cur_mgi, el_nd)
   cur_e_dens(cur_mgi) = el_nd
   cur_J(cur_mgi) = 0.D0   
  END IF ! test for the first iteration
  write(propmod_file,"(A, A12)") TRIM(outputfolder), '/propmod.dat'
  INQUIRE(FILE=propmod_file, EXIST=propmod_file_exists)

  ! write(*,*) 'update_grid: saved_grid = ', saved_grid, ' propmod_file_exists = ', propmod_file_exists
  IF(saved_grid == 1 .and. propmod_file_exists) THEN
  ! nothing
  ELSE
  ! if we calculate the condition only from electron density
   ! write(*,*) 'update_grid: enable_diffusion = ', enable_diffusion
   IF(enable_diffusion == 1) THEN
    CALL diffusion_approximation(cur_mgi, cur_e_dens(cur_mgi))
   END IF
  END IF
  ! write(*,*) 'update_grid: temp = ', model_grid(cur_mgi)%T
 END IF ! cur_n_assoccells > 0
END DO

! write(*,*) 'update_grid: my_rank = ', my_rank, ' the physical structure is computed'

#if mpi == 1

! DO cur_mgi = my_start, my_end
!  ! test_temp = model_grid(cur_mgi)%T
!  IF(model_grid(cur_mgi)%assoc_cells > 0) THEN
!   test_temp = cur_temp(cur_mgi)
!   write(*,*) 'update_grid: my_rank = ', my_rank, ' cur_temp = ', cur_temp(cur_mgi)
!   IF(test_temp == 0.D00) THEN
!    write(*,*) 'update_grid: before reduce'
!    write(*,*) 'update_grid: my_rank = ', my_rank, ' cur_mgi = ', cur_mgi
!    STOP 'update_grid: temp = 0 K'
!   END IF
!  END IF
! END DO

 IF(n_tasks > 1) THEN
  CALL MPI_ALLREDUCE(cur_j(1:n_modelgrid + add_mg), recv_j(:n_modelgrid + add_mg), n_modelgrid + add_mg, &
   MPI_DOUBLE, MPI_SUM, mpi_comm_world, ierr)
  CALL MPI_ALLREDUCE(cur_temp(1:n_modelgrid + add_mg), recv_temp(1:n_modelgrid + add_mg), n_modelgrid + add_mg, &
   MPI_DOUBLE, MPI_SUM, mpi_comm_world, ierr)
  IF(nlte /= 5) THEN
   CALL MPI_ALLREDUCE(cur_e_dens(1:n_modelgrid + add_mg), recv_e_dens(1:n_modelgrid + add_mg), n_modelgrid + add_mg, &
    MPI_DOUBLE, MPI_SUM, mpi_comm_world, ierr)
  END IF
  CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
  model_grid(:)%j = recv_j(:)
  IF(nlte /= 5) model_grid(:)%T = recv_temp(:)
  model_grid(:)%e_dens = recv_e_dens(:)
 ELSE IF(n_tasks == 1) THEN
  model_grid(1:n_modelgrid)%j = cur_j(1:n_modelgrid)
  IF(nlte /= 5) model_grid(1:n_modelgrid)%T = cur_temp(1:n_modelgrid)
  model_grid(1:n_modelgrid)%e_dens = cur_e_dens(1:n_modelgrid)
 END IF
#endif



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
