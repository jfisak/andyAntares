! connection propGrid and the modGrid for the basic case
!
! INPUT: NONE
! OUTPUT: NONE
!
SUBROUTINE connect_1D_basic()

USE types
USE MPI

IMPLICIT NONE

! loop variables
INTEGER                        :: cur_propcell, ind_J, best_index
DOUBLE PRECISION               :: delta, delta2
DOUBLE PRECISION               :: radius
DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

  

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek

INTEGER, DIMENSION(n_propgcells)                :: cur_model_index
INTEGER, DIMENSION(n_modelgrid + add_mg)        :: cur_n_assocmodg



#if mpi == 1
 N_single = (n_modelgrid)/n_tasks
 N_zbytek = n_modelgrid - n_tasks * N_single
 IF(my_rank <= N_zbytek - 1) THEN
  my_start = my_rank * (N_single + 1) + 1
  my_end = my_rank * (N_single + 1) + N_single
 ELSE IF(N_zbytek == 0) THEN
  my_start = my_rank * (N_single) + 1
  my_end = my_rank * (N_single) + N_single
 ELSE
  my_start = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + 1
  my_end = N_zbytek * (N_single + 1) + (my_rank - N_zbytek - 1) * N_single + N_single +1
 END IF

 IF(my_rank == n_tasks - 1) THEN
  my_end = n_modelgrid
 END IF
#else
 my_start = 1
 my_end = n_propgcells
#endif
write(*,*) 'connection_prop_model_grid: my_rank = ', my_rank, ' n_modelgrid = ', n_modelgrid
write(*,*) 'connection_prop_model_grid: N_single = ', N_single, ' N_zbytek = ', N_zbytek
write(*,*) 'connection_prop_model_grid: my_rank = ', my_rank, ' my_start = ', my_start, ' my_end = ', my_end


DO cur_propcell = my_start, my_end
 IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
  radius = SQRT((dyn_cell(cur_propcell)%corner(ind_x) + dyn_cell(cur_propcell)%width(ind_x)/2.D0)**2 + &
   (dyn_cell(cur_propcell)%corner(ind_y) + dyn_cell(cur_propcell)%width(ind_y)/2.D0)**2 + &
   (dyn_cell(cur_propcell)%corner(ind_z) + dyn_cell(cur_propcell)%width(ind_z)/2.D0)**2)
  ! IF ((radius .GT. R_star) .AND. (radius .LT. R_inf)) THEN
  ! Cells with radius larger than the stellar radius but smaller
  ! than the winds outer radius have an associated model grid cell.
  ! Find this model grid cell and add a pointer to the propatation
  ! grid. Finally record the number of asscociated prop. grid cells
  ! on the model grid
  delta = large_number
  best_index = 0
  DO ind_J = 1, n_modelgrid   
   delta2 = ABS(radius - model_grid(ind_J)%rwind)
   IF (delta2 .LT. delta) THEN
    delta = delta2 
    best_index = ind_J           
   END IF
  END DO
  IF(best_index > 0) THEN
   cur_model_index(cur_propcell) = best_index
   cur_n_assocmodg(best_index) = cur_n_assocmodg(best_index) + 1
  ELSE ! best_index <= 0
    cur_model_index(cur_propcell) = outerspace_index
   cur_n_assocmodg(outerspace_index) = cur_n_assocmodg(outerspace_index) + 1
  END IF ! best_index > 0
 END IF ! up_cell == 0
END DO ! a loop over propGrid cells

#if mpi == 1
 ! write(*,*) 'connection_prop_model_grid: ', SIZE(cur_model_index), SIZE(dyn_cell(:)%model_index), n_propgcells
 ! STOP 'connection_prop_model_grid: testing'
 CALL MPI_REDUCE(cur_model_index(:), dyn_cell(:)%model_index, n_propgcells, MPI_INTEGER, &
  & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_REDUCE(cur_n_assocmodg(:), model_grid(:)%assoc_cells, n_modelgrid, MPI_INTEGER, &
  & MPI_SUM, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_BCAST(dyn_cell(:)%model_index, n_propgcells, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
 CALL MPI_BCAST(model_grid(:)%assoc_cells, n_modelgrid + add_mg, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif 











END SUBROUTINE connect_1D_basic
