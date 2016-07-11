  SUBROUTINE setup_grid2() 
 
  ! Set up propagation grid cells

  USE types

  IMPLICIT NONE    

  INTEGER                        :: I, J, K, L, M
  INTEGER, PARAMETER             :: Nmax = 1000000000               
  DOUBLE PRECISION               :: r, delta, delta2           
  DOUBLE PRECISION, DIMENSION(3) :: cell_width2
  INTEGER                        :: max_n_dcell


  ! Number of propagation grid cells
  Ngrid = nx_cell * ny_cell * nz_cell

  IF (Ngrid .GT. Nmax) THEN
     PRINT*, 'ERROR: N > Nmax', Ngrid
     STOP 
  END IF

  ! Define length of dynamic arrays (cell and package)
  ALLOCATE (dyn_cell(2 * Ngrid))

  ! Size of the grid cells in x,y, and z direction (now they are with
  ! the same size i.e. regular gred)
  cell_width2(1) = 2.E0 * xmax / nx_cell
  cell_width2(2) = 2.E0 * ymax / ny_cell
  cell_width2(3) = 2.E0 * zmax / nz_cell

  L = 1
  DO I=1, nx_cell
   DO J=1, ny_cell
    DO K=1, nz_cell
     ! Index(number) of each cell in x,y, and z direction
!     dyn_cell(L)%indexc(1) = I
!     dyn_cell(L)%indexc(2) = J
!     dyn_cell(L)%indexc(3) = K
     ! Coordinates of the lower left corner of each cell
     dyn_cell(L)%corner(1)  = - xmax + (I - 1) * cell_width2(1)
     dyn_cell(L)%corner(2)  = - ymax + (J - 1) * cell_width2(2)     
     dyn_cell(L)%corner(3)  = - zmax + (K - 1) * cell_width2(3) 
     ! cell width
     dyn_cell(L)%width(1) = cell_width2(1)
     dyn_cell(L)%width(2) = cell_width2(2)
     dyn_cell(L)%width(3) = cell_width2(3)
     ! number of down cell is equal to zero
     dyn_cell(L)%down_cell = 0
     dyn_cell(L)%up_cell = 0
     write(15,*) dyn_cell(L)%corner, dyn_cell(L)%width
!       write(2,*) cell(L)%indexc, cell(L)%corner ! don't write if it is not necessary (computational very expensive)
     L = L + 1
    END DO
   END DO
  END DO
  CLOSE(15)
  ! the maximal number of cells is now equal to L
  max_n_dcell = Ngrid

  DO I = 1, Ngrid
   CALL create_dynamical_grid_cells(I, max_n_dcell)
   print*, 'max_n_dcell = ', max_n_dcell
  END DO
  OPEN(15,FILE='dyn_cells.dat')
  DO I = 1, max_n_dcell
   write(15,*) dyn_cell(I)%corner, dyn_cell(I)%width
  END DO


  
  ! Establish a connection between the propagation grid and the
  ! model grid. This depends on the model grid type (1D, 2D, 3D)
!  IF (model_type .EQ. 1) THEN
!    ! This is the algorithm needed for a 1D model grid
!    ! Define which model grid cell coresponds to the propagation grid cell
!    DO I = 1, Ngrid
!      ! Absolute radius of the propagation grid cell (midle of the cell)
!      r =SQRT( (cell(I)%corner(1) + cell_width/2.D0)**2 + &
!              (cell(I)%corner(2) + cell_width/2.D0)**2 + &
!              (cell(I)%corner(3) + cell_width/2.D0)**2)
!      !print*,I,r/R_star
!      IF ((r .GT. R_star) .AND. (r .LT. R_inf)) THEN
!        ! Cells with radius larger than the stellar radius but smaller
!        ! than the winds outer radius have an associated model grid cell.
!        ! Find this model grid cell and add a pointer to the propatation
!        ! grid. Finally record the number of asscociated prop. grid cells
!        ! on the model grid
!        delta = 1.D99
!        DO J = 1, n_modelgrid   
!           delta2 = ABS(r - model_grid(J)%rwind)
!           !print*,I,J,r/R_star,model_grid(J)%rwind/R_star,delta2/R_star,delta/R_star
!           IF (delta2 .LT. delta) THEN
!              delta = delta2 
!              M = J           
!           END IF
!        END DO
!        cell(I)%model_index = M     
!        model_grid(M)%assoc_cells = model_grid(M)%assoc_cells + 1
!      ELSE
!        ! Cells with radius smaller than the stellar radius or larger
!        ! than the winds outer radius have no associated model grid cell
!        ! Make them point to the dummy model grid cell
!        cell(I)%model_index = n_modelgrid + 1     
!      ENDIF
!      !print*, I,J,M
!    END DO
!  ENDIF
!
!
!  DO I = 1, n_modelgrid
!     print*, I, model_grid(I)%assoc_cells
!  END DO


END SUBROUTINE setup_grid2
