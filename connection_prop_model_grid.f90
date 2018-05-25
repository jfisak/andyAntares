  SUBROUTINE connection_prop_model_grid()

  USE types

  IMPLICIT NONE
  ! maximal distance between model and propagation grid
  ! MUST BE LATER CHANGED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DOUBLE PRECISION               :: basic_diagonal
  DOUBLE PRECISION               :: diagonal
  ! loop variables
  INTEGER                        :: I, J, M
  INTEGER                        :: max_n_dcell
  ! variables for calculating the shortest distance between
  ! propagation and model cell
  DOUBLE PRECISION               :: delta, delta2
  ! radial and vertical distance
  DOUBLE PRECISION               :: r, z, r0, z0
  ! volume of model cell
  DOUBLE PRECISION               :: volume, loc_volume
  INTEGER                        :: gridcell
  
  basic_diagonal = sqrt(basic_cell_width(1)**2 + basic_cell_width(2)**2 + &
                        basic_cell_width(3)**2)

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
      delta = 1.D99
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
   add_mg = 2
   !$OMP PARALLEL
   !$DEFAULT(private)
   !$OMP DO 
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
        !diagonal = sqrt(dyn_cell(I)%width(1)**2+dyn_cell(I)%width(2)**2+dyn_cell(I)%width(3)**2)
        diagonal = sqrt(dyn_cell(I)%width(1)**2+dyn_cell(I)%width(3)**2)/2.D0
        !diagonal = sqrt((sqrt(r**2 - z**2) - sqrt(r0**2 - z0**2))**2 + &
        ! (z - z0)**2)
        !IF( (delta > diagonal) .AND. (diagonal < basic_diagonal / 16.D0)) THEN
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
   !$OMP END DO
   !$OMP END PARALLEL
   write(99,*) 'number of propagation cells in vacuum: ', model_grid(n_modelgrid + add_mg)%assoc_cells
  ENDIF
! computing volume of model cells
DO gridcell = 1, n_modelgrid
 volume = 0.D0
 DO I = 1, max_n_dcell
  IF(dyn_cell(I)%up_cell == 0) THEN
   IF(dyn_cell(I)%model_index == gridcell) THEN
    loc_volume = dyn_cell(I)%width(1) * dyn_cell(I)%width(2) * dyn_cell(I)%width(3)
    volume = volume + loc_volume
   END IF
  END IF
 END DO
 !print*, 'connection_prop_model_grid: volume of the cell ', gridcell, ' is ', volume
 model_grid(gridcell)%volume = volume
END DO

!  print*, 'printing number of associated cells'
!  OPEN(UNIT=3,FILE='conneced_cells.dat')
!   DO I=1, n_modelgrid
!!    IF(dyn_cell(I)%model_index == n_modelgrid + add_mg) write(3,*) dyn_cell(I)%corner, dyn_cell(I)%width
!     write(3,*), I, model_grid(I)%assoc_cells
!   END DO
!  CLOSE(3)
!  DO I = 1, n_modelgrid + 1
!     print*, I, model_grid(I)%assoc_cells
!  END DO

  END SUBROUTINE connection_prop_model_grid
