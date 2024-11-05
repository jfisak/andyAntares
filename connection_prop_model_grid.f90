! connects propGrid and modGrid cells
! for each propGrid finds a corresponding modGrid cell
! the variable is saved into dyn_cell(:)%model_index (INT)
! and counted in the model_grid(:)%assoc_cells (INT)
! 1 ... n_modelgrid
! or
! photosphere_index
! outerspace_index
! vacuum_index
!
! INPUT: NONE
! OUTPUT: NONE
!
! the included function will be moved to special subroutines as well as the 3D case
SUBROUTINE connection_prop_model_grid()

USE MPI
USE types
USE constants
USE counters

IMPLICIT NONE
! maximal distance between model and propagation grid
! MUST BE LATER CHANGED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION               :: diagonal
! loop variables
INTEGER                        :: cur_propcell, J, best_index
! variables for calculating the shortest distance between
! propagation and model cell
DOUBLE PRECISION               :: delta, delta2
! radial and vertical distance
DOUBLE PRECISION               :: r, z, r0, phi, phi0

DOUBLE PRECISION               :: dist

DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

  

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek, N_tot_zbytek

INTEGER, DIMENSION(n_propgcells)                :: cur_model_index
INTEGER, DIMENSION(n_modelgrid + add_mg)        :: cur_n_assocmodg


cur_model_index(:) = 0
cur_n_assocmodg(:) = 0

! Establish a connection between the propagation grid and the
! model grid. This depends on the model grid type (1D, 2D, 3D)
!__________________________________________________________________________________________________
!__________________________________________________________________________________________________
! 1D case
!__________________________________________________________________________________________________
!__________________________________________________________________________________________________
IF (model_type .EQ. 1) THEN
 ! This is the algorithm needed for a 1D model grid
 ! Define which model grid cell coresponds to the propagation grid cell
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
     ! Absolute radius of the propagation grid cell (midle of the cell)
!      r = SQRT( (dyn_cell(I)%corner(1) + dyn_cell(I)%width(1)/2.D0)**2 + &
!       (dyn_cell(I)%corner(2) + dyn_cell(I)%width(2)/2.D0)**2 + &
!       (dyn_cell(I)%corner(3) + dyn_cell(I)%width(3)/2.D0)**2)
!      IF ((r .GT. R_star) .AND. (r .LT. R_inf)) THEN
!       ! Cells with radius larger than the stellar radius but smaller
!       ! than the winds outer radius have an associated model grid cell.
!       ! Find this model grid cell and add a pointer to the propatation
!       ! grid. Finally record the number of asscociated prop. grid cells
!       ! on the model grid
!       delta = large_number
!       DO J = 1, n_modelgrid   
!        delta2 = ABS(r - model_grid(J)%rwind)
!        IF (delta2 .LT. delta) THEN
!         delta = delta2 
!         M = J           
!        END IF
!       END DO
!       dyn_cell(I)%model_index = M     
!       model_grid(M)%assoc_cells = model_grid(M)%assoc_cells + 1
!      ELSE
!       ! Cells with radius smaller than the stellar radius or larger
!       ! than the winds outer radius have no associated model grid cell
!       ! Make them point to the dummy model grid cell
!       dyn_cell(I)%model_index = n_modelgrid + 1     
!       model_grid(n_modelgrid + 1)%assoc_cells = model_grid(n_modelgrid + 1)%assoc_cells + 1
!      END IF
!     END IF
!    END DO
     r = SQRT( (dyn_cell(cur_propcell)%corner(1) + dyn_cell(cur_propcell)%width(1)/2.D0)**2 + &
      (dyn_cell(cur_propcell)%corner(2) + dyn_cell(cur_propcell)%width(2)/2.D0)**2 + &
      (dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0)**2)
     ! IF ((r .GT. R_star) .AND. (r .LT. R_inf)) THEN
     ! Cells with radius larger than the stellar radius but smaller
     ! than the winds outer radius have an associated model grid cell.
     ! Find this model grid cell and add a pointer to the propatation
     ! grid. Finally record the number of asscociated prop. grid cells
     ! on the model grid
     delta = large_number
     best_index = 0
     DO J = 1, n_modelgrid   
      delta2 = ABS(r - model_grid(J)%rwind)
      IF (delta2 .LT. delta) THEN
       delta = delta2 
       best_index = J           
      END IF
     END DO
     IF(best_index > 0) THEN
      cur_model_index(cur_propcell) = best_index
      cur_n_assocmodg(best_index) = cur_n_assocmodg(best_index) + 1
     ELSE ! best_index <= 0
       cur_model_index(cur_propcell) = n_modelgrid + 1
      cur_n_assocmodg(n_modelgrid + 1) = cur_n_assocmodg(n_modelgrid + 1) + 1
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
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ! 2D case
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ELSE IF (model_type .EQ. 2) THEN
   SELECT CASE (inputmodel)
   ! basic 2D model
   CASE(1)
    CALL connect_2D_basic()
   ! PeKu model
   CASE(2)
    CALL connect_2D_peku()
    ! connect every single cell to its model cell
    DO cur_propcell = 1, n_propgcells
     r = SQRT((dyn_cell(cur_propcell)%corner(1) + dyn_cell(cur_propcell)%width(1)/2.D0)**2 + &
              (dyn_cell(cur_propcell)%corner(2) + dyn_cell(cur_propcell)%width(2)/2.D0)**2 + &
              (dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0)**2)
     z = dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0
     phi = acos(z/r)
     phi = abs(phi)
     IF(r < R_star .OR. r > R_inf) THEN
      dyn_cell(cur_propcell)%model_index = n_modelgrid
      CONTINUE
     END IF
     ! write(*,*) 'connection_prop_model_grid: r = ', r, ' z = ', z
     ! write(*,*) 'connection_prop_model_grid: phi = ', phi
      delta = 1.D99
      DO J = 1, n_modelgrid
       r0 = model_grid(J)%rwind
       phi0 = model_grid(J)%angle
       delta2 = sqrt(r**2.0+r0**2.0 - 2.0 * r * r0 * &
        (cos(phi)*cos(phi0) - sin(phi) * sin(phi0)))
       IF( delta2 < delta ) THEN
         delta = delta2
         best_index = J
        ! write(*,*) 'connection_prop_model_grid: cur_propcell = ', cur_propcell, ' / ', r/r0, phi/phi0
       END IF
       ! if the propagation cell is too far from the nearest model point
       ! we will associate this cell to the dummy cells
      END DO
       diagonal = sqrt(dyn_cell(cur_propcell)%width(1)**2+dyn_cell(cur_propcell)%width(3)**2)/2.D0
       IF((delta > diagonal) .AND. (dyn_cell(cur_propcell)%width(1) > basic_cell_width(1)/2.D0**6)) THEN
        dyn_cell(cur_propcell)%model_index = n_modelgrid + add_mg
        model_grid(n_modelgrid + add_mg)%assoc_cells = model_grid(n_modelgrid + add_mg)%assoc_cells + 1
       ELSE
        dyn_cell(cur_propcell)%model_index = best_index     
        IF(dyn_cell(cur_propcell)%model_index == 0) THEN
         write(*,*) 'connection_prop_model_grid: a cell ', cur_propcell, 'is not connected...'
         STOP
        END IF
        model_grid(best_index)%assoc_cells = model_grid(best_index)%assoc_cells + 1
       END IF
    END DO
    ! write(*,*) 'connection_prop_model_grid: my_rank = ', my_rank, ' my_start = ', my_start, &
    !  ' my_end = ', my_end
    ! STOP 'connection_prop_model_grid: testing'
   CASE DEFAULT
    STOP
   END SELECT
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ! 3D case
  !__________________________________________________________________________________________________
  !__________________________________________________________________________________________________
  ELSE IF (model_type == 3) THEN
   SELECT CASE(inputmodel)
   ! pseudo 3D testing model
   CASE(0)
    CALL connect_3D_pseudo()
   ! hydronico model
   CASE(1)
    CALL connect_3D_hydronico()
   CASE DEFAULT
    write(*,*) 'the choice inputmodel = ', inputmodel, ' is not known'
    STOP
   END SELECT
  END IF

END SUBROUTINE connection_prop_model_grid
