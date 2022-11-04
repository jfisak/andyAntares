SUBROUTINE connection_prop_model_grid()

USE types
USE counters

IMPLICIT NONE
! maximal distance between model and propagation grid
! MUST BE LATER CHANGED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DOUBLE PRECISION               :: diagonal
! loop variables
INTEGER                        :: I, J, M
INTEGER                        :: max_n_dcell
! variables for calculating the shortest distance between
! propagation and model cell
DOUBLE PRECISION               :: delta, delta2
! radial and vertical distance
DOUBLE PRECISION               :: r, z, r0, z0, phi, phi0
! volume of model cell
DOUBLE PRECISION               :: loc_volume
INTEGER                        :: gridcell
INTEGER                        :: my_n_cells
INTEGER                       :: N0, Nzbytek
INTEGER                       :: my_start, my_end, zb

DOUBLE PRECISION, DIMENSION(3) :: cur_center, cur_mpos
DOUBLE PRECISION               :: dist
INTEGER                        :: cur_mcell

DOUBLE PRECISION, PARAMETER     :: large_number=1.d90

DOUBLE PRECISION, DIMENSION(3)  :: cur_corner, cur_width
INTEGER                         :: cur_pgcell, cur_pgi
INTEGER                         :: n_virtpoints, cur_vpoint

DOUBLE PRECISION, DIMENSION(3)  :: cur_vpos
INTEGER                         :: up_cell, cur_lowcell
  
INTEGER, ALLOCATABLE            :: list_index(:)
INTEGER                         :: last_index, new_index
INTEGER                         :: n_bas_pcell

INTEGER                         :: ind_x, ind_y, ind_z
INTEGER                         :: cur_bpgi
DOUBLE PRECISION, DIMENSION(3)  :: width, cur_pos

INTEGER                         :: down_cell
INTEGER                         :: start_vp_index, end_vp_index
INTEGER                         :: cur_vp_nearest, cur_vp_index
INTEGER                         :: act_pgcell, cur_vp, cur_mgi


  max_n_dcell = SIZE(dyn_cell)
  ! Establish a connection between the propagation grid and the
  ! model grid. This depends on the model grid type (1D, 2D, 3D)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 1D model grid -- radial symetric
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (model_type .EQ. 1) THEN
   ! This is the algorithm needed for a 1D model grid
   ! Define which model grid cell coresponds to the propagation grid cell
   DO I = 1, max_n_dcell
    IF(dyn_cell(I)%up_cell == 0) THEN
     ! Absolute radius of the propagation grid cell (midle of the cell)
     r = SQRT( (dyn_cell(I)%corner(1) + dyn_cell(I)%width(1)/2.D0)**2 + &
      (dyn_cell(I)%corner(2) + dyn_cell(I)%width(2)/2.D0)**2 + &
      (dyn_cell(I)%corner(3) + dyn_cell(I)%width(3)/2.D0)**2)
     !print*,I,r/R_star
     IF ((r .GT. R_star) .AND. (r .LT. R_inf)) THEN
      ! Cells with radius larger than the stellar radius but smaller
      ! than the winds outer radius have an associated model grid cell.
      ! Find this model grid cell and add a pointer to the propatation
      ! grid. Finally record the number of asscociated prop. grid cells
      ! on the model grid
      delta = large_number
      DO J = 1, n_modelgrid   
       delta2 = ABS(r - model_grid(J)%rwind)
       !print*,I,J,r/R_star,model_grid(J)%rwind/R_star,delta2/R_star,delta/R_star
       IF (delta2 .LT. delta) THEN
        delta = delta2 
        M = J           
       END IF
      END DO
      dyn_cell(I)%model_index = M     
      model_grid(M)%assoc_cells = model_grid(M)%assoc_cells + 1
     ELSE
      ! Cells with radius smaller than the stellar radius or larger
      ! than the winds outer radius have no associated model grid cell
      ! Make them point to the dummy model grid cell
      dyn_cell(I)%model_index = n_modelgrid + 1     
      model_grid(n_modelgrid + 1)%assoc_cells = model_grid(n_modelgrid + 1)%assoc_cells + 1
     END IF
    END IF
     !print*, I,J,M
   END DO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 2D model grid -- Petr Kurfurst's model
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ELSE IF (model_type .EQ. 2) THEN
   SELECT CASE (inputmodel)
   ! PeKu disk model
   CASE(1)
   add_mg = 2
   DO I = 1, max_n_dcell
    ! IF(mod(I,10000) .EQ. 0) print*, 'associating propagation grid', I, REAL(I)/REAL(max_n_dcell) * 1.E2, ' % completed'
    IF(dyn_cell(I)%up_cell == 0) THEN
      ! Absolute radius of the propagation grid cell (midle of the cell)
      r = SQRT((dyn_cell(I)%corner(1) + dyn_cell(I)%width(1)/2.D0)**2 + &
               (dyn_cell(I)%corner(2) + dyn_cell(I)%width(2)/2.D0)**2 + &
               (dyn_cell(I)%corner(3) + dyn_cell(I)%width(3)/2.D0)**2)
      r0 = SQRT(dyn_cell(I)%corner(1)**2 + dyn_cell(I)%corner(2)**2 + &
               dyn_cell(I)%corner(3)**2)
      z = dyn_cell(I)%corner(3) + dyn_cell(I)%width(3)/2.D0
      z0 = dyn_cell(I)%corner(3)
      !print*,I,r/R_star
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
          M = J
        END IF
        ! if the propagation cell is too far from the nearest model point
        ! we will associate this cell to the dummy cells
       END DO
        diagonal = sqrt(dyn_cell(I)%width(1)**2+dyn_cell(I)%width(3)**2)/2.D0
        IF((delta > diagonal) .AND. (dyn_cell(I)%width(1) > basic_cell_width(1)/2.D0**6)) THEN
         dyn_cell(I)%model_index = n_modelgrid + add_mg
         model_grid(n_modelgrid + add_mg)%assoc_cells = model_grid(n_modelgrid + add_mg)%assoc_cells + 1
         !print*, 'model grid n + 2 = ', model_grid(n_modelgrid + 2)%assoc_cells
        ELSE
         dyn_cell(I)%model_index = M     
         model_grid(M)%assoc_cells = model_grid(M)%assoc_cells + 1
        END IF
      ELSE
       ! Cells with radius smaller than the stellar radius or larger
       ! than the winds outer radius have no associated model grid cell
       ! Make them point to the dummy model grid cell
       dyn_cell(I)%model_index = n_modelgrid + 1     
       !model_grid(n_modelgrid + 1)%assoc_cells = model_grid(n_modelgrid + 1)%assoc_cells + 1
      END IF
      !print*, I,J,M
    END IF
   END DO
   write(99,*) 'number of propagation cells in vacuum: ', model_grid(n_modelgrid + add_mg)%assoc_cells
   ! supernova model
   CASE(2)
    ! connect every single cell to its model cell
    DO I = my_start, my_end
     r = SQRT((dyn_cell(I)%corner(1) + dyn_cell(I)%width(1)/2.D0)**2 + &
              (dyn_cell(I)%corner(2) + dyn_cell(I)%width(2)/2.D0)**2 + &
              (dyn_cell(I)%corner(3) + dyn_cell(I)%width(3)/2.D0)**2)
     z = dyn_cell(I)%corner(3) + dyn_cell(I)%width(3)/2.D0
     phi = acos(z/r)
     phi = abs(phi)
     IF(r < R_star .OR. r > R_inf) THEN
      dyn_cell(I)%model_index = n_modelgrid
      CONTINUE
     END IF
     write(*,*) 'connection_prop_model_grid: r = ', r, ' z = ', z
     write(*,*) 'connection_prop_model_grid: phi = ', phi
      delta = 1.D99
      DO J = 1, n_modelgrid
       r0 = model_grid(J)%rwind
       phi0 = model_grid(J)%angle
       delta2 = sqrt(r**2.0+r0**2.0 - 2.0 * r * r0 * &
        (cos(phi)*cos(phi0) - sin(phi) * sin(phi0)))
       IF( delta2 < delta ) THEN
         delta = delta2
         M = J
       END IF
       ! if the propagation cell is too far from the nearest model point
       ! we will associate this cell to the dummy cells
      END DO
       diagonal = sqrt(dyn_cell(I)%width(1)**2+dyn_cell(I)%width(3)**2)/2.D0
       IF((delta > diagonal) .AND. (dyn_cell(I)%width(1) > basic_cell_width(1)/2.D0**6)) THEN
        dyn_cell(I)%model_index = n_modelgrid + add_mg
        model_grid(n_modelgrid + add_mg)%assoc_cells = model_grid(n_modelgrid + add_mg)%assoc_cells + 1
       ELSE
        dyn_cell(I)%model_index = M     
        model_grid(M)%assoc_cells = model_grid(M)%assoc_cells + 1
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
    DO I = 1, max_n_dcell
     dyn_cell(I)%model_index = I
     model_grid(I)%assoc_cells = model_grid(I)%assoc_cells + 1
    END DO
   CASE(1)
    ! create an array with saved indexes
    n_bas_pcell = nx_cell * ny_cell * nz_cell
    ALLOCATE(list_index(n_bas_pcell))
    ! initial setup
    DO I = 1, n_bas_pcell
     list_index(I) = 0
    END DO
    ! calculating of indeces
    last_index = 0
    n_virtpoints = SIZE(virtual_point)
    DO I = 1, n_virtpoints
     new_index = virtual_point(I)%ind_pcell
     IF(new_index /= last_index) THEN
      list_index(new_index) = I
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
!        write(*,*) 'connection_prop_model_grid: cur_pgcell = ', cur_pgcell, ' modindex = ', virtual_point(cur_vp_nearest)%ind_mcell

       END IF
      END IF
     END IF ! up_cell == 0
    END DO 
   CASE DEFAULT
    write(*,*) 'the choice inputmodel = ', inputmodel, ' is not known'
    STOP
   END SELECT
  END IF
! computing volume of model cells
DO I = 1, max_n_dcell
 IF(dyn_cell(I)%up_cell == 0) THEN
  gridcell = dyn_cell(I)%model_index
   loc_volume = dyn_cell(I)%width(1) * dyn_cell(I)%width(2) * dyn_cell(I)%width(3)
   model_grid(gridcell)%volume = model_grid(gridcell)%volume + loc_volume
  ! counters
  IF(dyn_cell(I)%model_index /= 0) THEN
   count_pg_mcell = count_pg_mcell + 1
  ELSE IF(dyn_cell(I)%model_index == n_modelgrid + 3) THEN
   count_pg_vacuum = count_pg_vacuum + 1
  END IF
 END IF
END DO

write(99,*) 'propGrid statistics:'
write(99,*) 'count_pg_mcell = ', count_pg_mcell
write(99,*) 'count_pg_vacuum = ', count_pg_vacuum

END SUBROUTINE connection_prop_model_grid
