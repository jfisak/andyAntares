  SUBROUTINE setup_grid2() 
 
  ! Set up propagation grid cells using dynamic cells

  USE types

  IMPLICIT NONE    

  ! loop variables
  INTEGER                                :: I, J, K, L, M
  INTEGER, PARAMETER                     :: Nmax = 1000000000               
  ! variables describing dynamic cells
  DOUBLE PRECISION, DIMENSION(3)         :: cell_width2
  INTEGER                                :: max_n_dcell, N_dyn_grid
  TYPE(dyn_grid_cell), ALLOCATABLE       :: pom2(:)


  ! Number of propagation grid cells
  Ngrid = nx_cell * ny_cell * nz_cell

  IF (Ngrid .GT. Nmax) THEN
     PRINT*, 'ERROR: N > Nmax', Ngrid
     STOP 
  END IF

  ! Define length of dynamic arrays (cell and package)
  ALLOCATE (dyn_cell(2 * Ngrid))

  ! Size of the basic grid cells in x,y, and z direction (now they are with
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
     !write(2,*) dyn_cell(L)%indexc, cell(L)%corner ! don't write if it is not necessary (computational very expensive)
     L = L + 1
    END DO
   END DO
  END DO
  CLOSE(15)
  ! the maximal number of cells is now equal to L
  max_n_dcell = Ngrid

 IF(dyngrid /= 0) THEN
  DO I = 1, Ngrid
   CALL create_dynamical_grid_cells(I, max_n_dcell)
   print*, 'max_n_dcell = ', max_n_dcell
  END DO
 END IF
  OPEN(15,FILE='dyn_cells.dat')
   DO I = 1, max_n_dcell
    write(15,*) dyn_cell(I)%corner, dyn_cell(I)%width
   END DO
  CLOSE(15)

   ! we will resize the field dyn_cell
   ! because we do not want empty cells
   ! inside
   ALLOCATE(pom2(max_n_dcell))
   do I = 1, max_n_dcell
    pom2(I) = dyn_cell(I)
   end do
   DEALLOCATE(dyn_cell)
   ALLOCATE(dyn_cell(max_n_dcell))
   do I = 1, max_n_dcell
    dyn_cell(I) = pom2(I)
   end do
   DEALLOCATE(pom2)
   N_dyn_grid = max_n_dcell

END SUBROUTINE setup_grid2
