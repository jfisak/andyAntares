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

INTEGER                         :: last_index

INTEGER                         :: cur_neighbour

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek
INTEGER                                 :: cur_mgi

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

add_mg = 2
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
     END IF
     ! if the propagation cell is too far from the nearest model point
     ! we will associate this cell to the dummy cells
    END DO ! loop over modGrid cells

     diagonal = sqrt(dyn_cell(cur_propcell)%width(1)**2+dyn_cell(cur_propcell)%width(3)**2)/2.D0
     IF((delta > diagonal) .AND. (dyn_cell(cur_propcell)%width(1) > basic_cell_width(1)/2.D0**6)) THEN
      dyn_cell(cur_propcell)%model_index = n_modelgrid + add_mg
      model_grid(n_modelgrid + add_mg)%assoc_cells = model_grid(n_modelgrid + add_mg)%assoc_cells + 1
     ELSE
      dyn_cell(cur_propcell)%model_index = best_index     
      model_grid(best_index)%assoc_cells = model_grid(best_index)%assoc_cells + 1
     END IF
   ELSE
    ! Cells with radius smaller than the stellar radius or larger
    ! than the winds outer radius have no associated model grid cell
    ! Make them point to the dummy model grid cell
    dyn_cell(cur_propcell)%model_index = n_modelgrid + 1     
    !model_grid(n_modelgrid + 1)%assoc_cells = model_grid(n_modelgrid + 1)%assoc_cells + 1
   END IF
   ! write(*,*) 'connect_2D_basic: cur_propcell = ', cur_propcell
   ! write(*,*) 'connect_2D_basic: modGrid index = ', dyn_cell(cur_propcell)%model_index
 END IF
END DO ! loop over propGrid cells
write(99,*) 'number of propagation cells in vacuum: ', model_grid(n_modelgrid + add_mg)%assoc_cells
END SUBROUTINE connect_2D_basic
