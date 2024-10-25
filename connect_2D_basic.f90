! a sbr connecting propGrid with the modGrid
! for a 2D model with defined
!
! radius and a lateral coordinate
!
! input: NONE
! output: NONE
!
SUBROUTINE connect_2D_basic()

USE types
USE constants

IMPLICIT NONE

DOUBLE PRECISION               :: diagonal
! loop variables
INTEGER                        :: cur_propcell, best_index
! variables for calculating the shortest distance between
! propagation and model cell
DOUBLE PRECISION               :: delta, delta2
! radial and vertical distance
DOUBLE PRECISION               :: pgi_radius, pgi_theta, pgi_radius0, pgi_theta0
DOUBLE PRECISION                                :: mgi_radius, mgi_theta
DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek
INTEGER                                 :: cur_mgi
DOUBLE PRECISION                :: tot_delta = 0.D0
INTEGER                         :: n_adjonced = 0

#if mpi == 1
   N_single = n_propgcells/n_tasks
   N_zbytek = n_propgcells - n_tasks * N_single
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
#else
   my_start = 1
   my_end = n_propgcells
#endif
my_start = 1
my_end = n_propgcells
DO cur_propcell = my_start, my_end
 ! IF(mod(cur_propcell,10000) .EQ. 0) print*, 'associating propagation grid', cur_propcell, REAL(cur_propcell)/REAL(max_n_dcell) * 1.E2, ' % completed'
 IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
   ! Absolute radius of the propagation grid cell (midle of the cell)
   pgi_radius = SQRT((dyn_cell(cur_propcell)%corner(1) + dyn_cell(cur_propcell)%width(1)/2.D0)**2 + &
            (dyn_cell(cur_propcell)%corner(2) + dyn_cell(cur_propcell)%width(2)/2.D0)**2 + &
            (dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0)**2)
   pgi_radius0 = SQRT(dyn_cell(cur_propcell)%corner(1)**2 + dyn_cell(cur_propcell)%corner(2)**2 + &
            dyn_cell(cur_propcell)%corner(3)**2)
   pgi_theta = acos((dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0)/pgi_radius)
   pgi_theta0 = acos(dyn_cell(cur_propcell)%corner(3)/pgi_radius0)



   IF ((pgi_radius .GT. R_star) .AND. (pgi_radius .LT. R_inf)) THEN
    ! Cells with radius larger than the stellar radius but smaller
    ! than the winds outer radius have an associated model grid cell.
    ! Find this model grid cell and add a pointer to the propatation
    ! grid. Finally record the number of asscociated prop. grid cells
    ! on the model grid
    delta = large_number
    DO cur_mgi = 1, n_modelgrid
     mgi_radius = model_grid(cur_mgi)%rwind
     mgi_theta = model_grid(cur_mgi)%angle
     ! calculation of distance
     delta2 = sqrt(mgi_radius**2 + pgi_radius**2 - 2.0 * mgi_radius * pgi_radius * cos(pgi_theta - mgi_theta))
     IF( delta2 < delta ) THEN
       delta = delta2
       best_index = cur_mgi
       n_adjonced = n_adjonced + 1
       tot_delta = tot_delta + delta2
     END IF
     ! if the propagation cell is too far from the nearest model point
     ! we will associate this cell to the dummy cells
    END DO ! loop over modGrid cells

     diagonal = sqrt(dyn_cell(cur_propcell)%width(1)**2+dyn_cell(cur_propcell)%width(3)**2)/2.D0
     IF((delta > diagonal) .AND. (dyn_cell(cur_propcell)%width(1) > basic_cell_width(1)/2.D0**6)) THEN
      dyn_cell(cur_propcell)%model_index = vacuum_index
      model_grid(vacuum_index)%assoc_cells = model_grid(vacuum_index)%assoc_cells + 1
     ELSE
      dyn_cell(cur_propcell)%model_index = best_index     
      model_grid(best_index)%assoc_cells = model_grid(best_index)%assoc_cells + 1
     END IF
   ELSE IF(pgi_radius < R_star) THEN
    ! Cells with radius smaller than the stellar radius or larger
    ! than the winds outer radius have no associated model grid cell
    ! Make them point to the dummy model grid cell
    dyn_cell(cur_propcell)%model_index = photosphere_index     
    model_grid(photosphere_index)%assoc_cells = model_grid(photosphere_index)%assoc_cells + 1
   ELSE IF(pgi_radius > R_inf) THEN
    dyn_cell(cur_propcell)%model_index = outerspace_index     
    model_grid(outerspace_index)%assoc_cells = model_grid(outerspace_index)%assoc_cells + 1
   END IF
   ! write(*,*) 'connect_2D_basic: cur_propcell = ', cur_propcell
   ! write(*,*) 'connect_2D_basic: modGrid index = ', dyn_cell(cur_propcell)%model_index
 END IF
END DO ! loop over propGrid cells
write(99,*) 'number of propagation cells in vacuum: ', model_grid(vacuum_index)%assoc_cells

DO cur_propcell = 1, n_propgcells
 IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
  IF(dyn_cell(cur_propcell)%model_index == 0) THEN
   write(*,*) 'connect_2D_basic: error, cur_propcell = ', cur_propcell
   STOP 'connect_2D_basic: model_index = 0'
  END IF
 END IF
END DO



END SUBROUTINE connect_2D_basic
