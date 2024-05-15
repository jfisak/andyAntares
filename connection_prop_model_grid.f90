SUBROUTINE connection_prop_model_grid()

USE types
USE constants
USE counters

IMPLICIT NONE
! maximal distance between model and propagation grid
! MUST BE LATER CHANGED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION               :: diagonal
! loop variables
INTEGER                        :: cur_propcell, J, best_index
INTEGER                        :: max_n_dcell
! variables for calculating the shortest distance between
! propagation and model cell
DOUBLE PRECISION               :: delta, delta2
! radial and vertical distance
DOUBLE PRECISION               :: r, z, r0, z0, phi, phi0
! volume of model cell
DOUBLE PRECISION               :: loc_volume
INTEGER                        :: gridcell

DOUBLE PRECISION, DIMENSION(3) :: cur_center
DOUBLE PRECISION               :: dist
INTEGER                        :: cur_mcell

DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

DOUBLE PRECISION, DIMENSION(3)  :: cur_corner, cur_width
INTEGER                         :: cur_pgcell, cur_pgi
INTEGER                         :: n_virtpoints

INTEGER                         :: up_cell, cur_lowcell
  
INTEGER, ALLOCATABLE            :: list_index(:)
INTEGER                         :: last_index, new_index
INTEGER                         :: n_bas_pcell

INTEGER                         :: cur_bpgi
DOUBLE PRECISION, DIMENSION(3)  :: cur_pos

INTEGER                         :: down_cell
INTEGER                         :: start_vp_index, end_vp_index
INTEGER                         :: cur_vp_nearest, cur_vp_index
INTEGER                         :: act_pgcell, cur_vp, cur_mgi

INTEGER                         :: count_vacuum, count_out, count_in, count_ok
INTEGER                         :: cur_neighbour, cur_n_mgi

! parallelization
INTEGER                         :: my_start, my_end
INTEGER                         :: N_single, N_zbytek

INTEGER, DIMENSION(n_propgcells)                :: cur_model_index
INTEGER, DIMENSION(n_modelgrid + add_mg)        :: cur_n_assocmodg

LOGICAL                                         :: vacuum_found

count_ok = 0
count_in = 0
count_out = 0
count_vacuum = 0

cur_model_index(:) = 0
cur_n_assocmodg(:) = 0

  max_n_dcell = SIZE(dyn_cell)
  ! Establish a connection between the propagation grid and the
  ! model grid. This depends on the model grid type (1D, 2D, 3D)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 1D model grid -- radial symetric
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (model_type .EQ. 1) THEN
   ! This is the algorithm needed for a 1D model grid
   ! Define which model grid cell coresponds to the propagation grid cell
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
write(*,*) 'update_grid: N_single = ', N_single, ' N_zbytek = ', N_zbytek
write(*,*) 'update_grid: my_rank = ', my_rank, ' my_start = ', my_start, ' my_end = ', my_end


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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2D model grid -- Petr Kurfurst's model
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE IF (model_type .EQ. 2) THEN
   SELECT CASE (inputmodel)
   ! PeKu disk model
   CASE(1)
   add_mg = 2
   DO cur_propcell = my_start, my_end
    ! IF(mod(cur_propcell,10000) .EQ. 0) print*, 'associating propagation grid', cur_propcell, REAL(cur_propcell)/REAL(max_n_dcell) * 1.E2, ' % completed'
    IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
      ! Absolute radius of the propagation grid cell (midle of the cell)
      r = SQRT((dyn_cell(cur_propcell)%corner(1) + dyn_cell(cur_propcell)%width(1)/2.D0)**2 + &
               (dyn_cell(cur_propcell)%corner(2) + dyn_cell(cur_propcell)%width(2)/2.D0)**2 + &
               (dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0)**2)
      r0 = SQRT(dyn_cell(cur_propcell)%corner(1)**2 + dyn_cell(cur_propcell)%corner(2)**2 + &
               dyn_cell(cur_propcell)%corner(3)**2)
      z = dyn_cell(cur_propcell)%corner(3) + dyn_cell(cur_propcell)%width(3)/2.D0
      z0 = dyn_cell(cur_propcell)%corner(3)
      !print*,cur_propcell,r/R_star
      IF ((r .GT. R_star) .AND. (r .LT. R_inf)) THEN
       ! Cells with radius larger than the stellar radius but smaller
       ! than the winds outer radius have an associated model grid cell.
       ! Find this model grid cell and add a pointer to the propatation
       ! grid. Finally record the number of asscociated prop. grid cells
       ! on the model grid
       delta = 1.D99
       DO J = 1, n_modelgrid
        delta2 = sqrt((sqrt(r**2 - z**2) - &
                (sqrt(model_grid(J)%rwind**2 - model_grid(J)%zwind**2)))**2 + &
                (model_grid(J)%zwind - z)**2)
        IF( delta2 < delta ) THEN
          delta = delta2
          best_index = J
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
         model_grid(best_index)%assoc_cells = model_grid(best_index)%assoc_cells + 1
        END IF
      ELSE
       ! Cells with radius smaller than the stellar radius or larger
       ! than the winds outer radius have no associated model grid cell
       ! Make them point to the dummy model grid cell
       dyn_cell(cur_propcell)%model_index = n_modelgrid + 1     
       !model_grid(n_modelgrid + 1)%assoc_cells = model_grid(n_modelgrid + 1)%assoc_cells + 1
      END IF
    END IF
   END DO
   write(99,*) 'number of propagation cells in vacuum: ', model_grid(n_modelgrid + add_mg)%assoc_cells
   ! PeKu model
   CASE(2)
    CALL connect_2D_peku()
    ! connect every single cell to its model cell
    DO cur_propcell = 1, max_n_dcell
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
  ELSE IF (model_type == 3) THEN
   SELECT CASE(inputmodel)
   ! pseudo 3D testing model
   CASE(0)
    DO cur_propcell = 1, max_n_dcell
     dyn_cell(cur_propcell)%model_index = cur_propcell
     model_grid(cur_propcell)%assoc_cells = model_grid(cur_propcell)%assoc_cells + 1
    END DO
   CASE(1)
    ! create an array with saved indexes
    ! write(*,*) 'connection_prop_model_grid: inputmodel = ', inputmodel, ' dyngrid = ', dyngrid
    IF(dyngrid > 0) THEN
     n_bas_pcell = nx_cell * ny_cell * nz_cell
     ALLOCATE(list_index(n_bas_pcell))
     ! initial setup
     DO cur_propcell = 1, n_bas_pcell
      list_index(cur_propcell) = 0
     END DO
     ! calculating of indeces
     last_index = 0
     n_virtpoints = SIZE(virtual_point)
     DO cur_propcell = 1, n_virtpoints
      new_index = virtual_point(cur_propcell)%ind_pcell
      IF(new_index /= last_index) THEN
       list_index(new_index) = cur_propcell
       last_index = new_index
      END IF
     END DO

     ! we calculate associated cells for the rest of propGrid cells
     DO cur_pgcell = 1, max_n_dcell
      up_cell = dyn_cell(cur_pgcell)%up_cell
      IF(up_cell == 0) THEN
       cur_corner = dyn_cell(cur_pgcell)%corner
       cur_width = dyn_cell(cur_pgcell)%width
       cur_center = cur_corner + cur_width/2.0
       cur_lowcell = dyn_cell(cur_pgcell)%down_cell
       r0 = sqrt(cur_center(1)**2.0 + cur_center(2)**2.0 + &
        cur_center(3)**2.0)
       ! add_mg = 1 r < R_star
       ! add_mg = 2 r > R_inf
       ! add_mg = 3 r > R_star && r < R_inf, vacuum cell
       IF(r0 < R_star) THEN
        dyn_cell(cur_pgcell)%model_index = n_modelgrid + 1
       ELSE IF(r0 > R_inf) THEN
        dyn_cell(cur_pgcell)%model_index = n_modelgrid + 2
       ELSE 
        ! we connect a modgrid from the current basic cell
        act_pgcell = cur_pgcell
        DO
         down_cell = dyn_cell(act_pgcell)%down_cell
         IF(down_cell == 0) EXIT
         act_pgcell = down_cell
        END DO
        cur_bpgi = act_pgcell

        start_vp_index = list_index(cur_bpgi)
        

        IF(start_vp_index == 0) THEN
         dyn_cell(cur_pgcell)%model_index = n_modelgrid + 3
        ELSE
         ! we must find an end index
         cur_vp_index = start_vp_index
         DO
          cur_vp_index = cur_vp_index + 1
          if(cur_vp_index == n_virtpoints) then
           end_vp_index = cur_vp_index
           EXIT
          end if
          if(virtual_point(cur_vp_index)%ind_pcell /= cur_bpgi) then
           end_vp_index = cur_vp_index - 1
           EXIT
          end if
         END DO

         delta = large_number
         cur_vp_nearest = 0
         DO cur_vp = start_vp_index, end_vp_index
          cur_pos = virtual_point(cur_vp)%pos
          dist = sqrt((cur_center(1) - cur_pos(1))**2.0 +&
           (cur_center(2) - cur_pos(2))**2.0 +&
           (cur_center(3) - cur_pos(3))**2.0)
          if(dist< delta) then
           delta = dist
           cur_vp_nearest = cur_vp
          end if
         END DO
         cur_mgi = virtual_point(cur_vp_nearest)%ind_mcell
         dyn_cell(cur_pgcell)%model_index = cur_mgi
         model_grid(cur_mgi)%assoc_cells = model_grid(cur_mgi)%assoc_cells + 1
!         write(*,*) 'connection_prop_model_grid: cur_pgcell = ', cur_pgcell, ' modindex = ', virtual_point(cur_vp_nearest)%ind_mcell

        END IF
       END IF
      END IF ! up_cell == 0
     END DO 
    ELSE IF(dyngrid == 0) THEN
     DO cur_propcell = 1, n_modelgrid
      cur_pos = model_grid(cur_propcell)%vec_pos
      CALL find_dyn_cell1(cur_pos, cur_pgi)
      IF(dyn_cell(cur_pgi)%model_index == 0) THEN
       dyn_cell(cur_pgi)%model_index = cur_propcell
       model_grid(cur_propcell)%assoc_cells = model_grid(cur_propcell)%assoc_cells + 1
       count_ok = count_ok + 1
      END IF
      DO J = 1,6
       cur_neighbour = dyn_cell(cur_pgi)%neighbor(J)
       IF(cur_neighbour > 0) then
        cur_n_mgi = dyn_cell(cur_neighbour)%model_index
        IF(cur_n_mgi == 0) THEN
         dyn_cell(cur_neighbour)%model_index = cur_propcell
         model_grid(cur_propcell)%assoc_cells = model_grid(cur_propcell)%assoc_cells + 1
         count_ok = count_ok + 1
        END IF
       END IF
      END DO
     END DO
     ! 
     vacuum_found = .true.
     DO WHILE(vacuum_found)
      vacuum_found = .false.
      DO cur_pgi = 1, max_n_dcell
       cur_mcell = dyn_cell(cur_pgi)%model_index
       if(cur_mcell == 0) then
        vacuum_found = .true.
       else
        DO J = 1,6
         IF(J == 5) cycle
         cur_neighbour = dyn_cell(cur_pgi)%neighbor(J)
         if(cur_neighbour > 0) then
          cur_n_mgi = dyn_cell(cur_neighbour)%model_index
          if(cur_n_mgi == 0) then
           dyn_cell(cur_neighbour)%model_index = cur_mcell
           model_grid(cur_mcell)%assoc_cells = model_grid(cur_mcell)%assoc_cells + 1
           ! write(*,*) 'connection_prop_model_grid: cur_pgcell = ', cur_pgcell, ' cur_mcell = ', cur_mcell
          end if ! cur_n_mgi == 0
         end if ! cur_neighbour > 0
        END DO 
       end if ! cur_mcell == 0
      END DO
     END DO
     ! other cells
     DO cur_pgcell = 1, max_n_dcell
      cur_corner = dyn_cell(cur_pgcell)%corner/R_star
      cur_width = dyn_cell(cur_pgcell)%width/R_star
      cur_center = cur_corner + cur_width/2.0
      r0 = sqrt(cur_center(1)**2.0 + cur_center(2)**2.0 + &
       cur_center(3)**2.0)
      IF(r0 < 1.0) THEN
       dyn_cell(cur_pgcell)%model_index = n_modelgrid + 1
       count_in = count_in + 1
      ELSE IF(r0 > R_inf/R_star) THEN
       dyn_cell(cur_pgcell)%model_index = n_modelgrid + 2
       count_out = count_out + 1
      END IF
      IF(dyn_cell(cur_pgcell)%model_index == 0) THEN
       dyn_cell(cur_pgcell)%model_index = n_modelgrid + 3
       count_vacuum = count_vacuum + 1
       ! delta2 = large_number
       ! DO J = 1, n_modelgrid
       !  mod_pos = model_grid(J)%vec_pos
       !  delta = sqrt((cur_center(1)-mod_pos(1))**2.0 + (cur_center(2)-mod_pos(2))**2.0 + (cur_center(3)-mod_pos(3))**2.0)
       !  if(delta < delta2) THEN
       !   delta2 = delta
       !   best_index = J
       !  end if
       !  if(mod(cur_propcell,1000)==0) write(*,*) 'connection_prop_model_grid: working on propGrid cell ', cur_propcell, ' from ', max_n_dcell
       ! END DO ! loop over all modGrid cells
       ! dyn_cell(cur_pgcell)%model_index = M
      END IF ! if model_index == 0
     END DO ! loop over all propGrid cells
     ! write(*,*) 'connection_prop_model_grid: in = ', count_in, ' out = ', count_out, ' vacuum = ', count_vacuum, &
     !  ' count_ok = ', count_ok
    END IF ! dyncell > 0
   CASE DEFAULT
    write(*,*) 'the choice inputmodel = ', inputmodel, ' is not known'
    STOP
   END SELECT
  END IF
! computing volume of model cells
DO cur_propcell = 1, max_n_dcell
 IF(dyn_cell(cur_propcell)%up_cell == 0) THEN
  gridcell = dyn_cell(cur_propcell)%model_index
   loc_volume = dyn_cell(cur_propcell)%width(1) * dyn_cell(cur_propcell)%width(2) * dyn_cell(cur_propcell)%width(3)
   model_grid(gridcell)%volume = model_grid(gridcell)%volume + loc_volume
  ! counters
  IF(dyn_cell(cur_propcell)%model_index /= 0) THEN
   count_pg_mcell = count_pg_mcell + 1
  ELSE IF(dyn_cell(cur_propcell)%model_index == n_modelgrid + 3) THEN
   count_pg_vacuum = count_pg_vacuum + 1
  END IF
 END IF
END DO

write(99,*) 'propGrid statistics:'
write(99,*) 'count_pg_mcell = ', count_pg_mcell
write(99,*) 'count_pg_vacuum = ', count_pg_vacuum

END SUBROUTINE connection_prop_model_grid
