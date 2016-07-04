  SUBROUTINE setup_basic_grid() 
 
  ! Set up propagation grid cells

  USE types

  IMPLICIT NONE    

  INTEGER             :: I, J, K, L, M
  INTEGER, PARAMETER  :: Nmax = 1000000000               
  DOUBLE PRECISION    :: r, delta, delta2           


  ! Number of propagation grid cells
  Ngrid = nx_cell * ny_cell * nz_cell

  IF (Ngrid .GT. Nmax) THEN
     PRINT*, 'ERROR: N > Nmax', Ngrid
     STOP 
  END IF

  ! Define length of dynamic arrays (cell and package)
  ALLOCATE (cell(Ngrid))

  ! Size of the grid cells in x,y, and z direction (now they are with
  ! the same size i.e. regular grid)
  cell_width = 2.D0 * xmax / nx_cell
  print*, 'parameters of grid: '
  print*, xmax, R_inf, cell_width
  print*, xmax/R_star, R_inf / R_star, cell_width / R_star
  print*, (-xmax + nx_cell * cell_width) / R_star

  L = 1
  DO I=1, nx_cell
     DO J=1, ny_cell
        DO K=1, nz_cell
           ! Index(number) of each cell in x,y, and z direction
           cell(L)%indexc(1) = I
           cell(L)%indexc(2) = J
           cell(L)%indexc(3) = K
!           ! Size of each grid cells in x,y, and z direction
!           cell(L)%deltax = deltax
!           cell(L)%deltay = deltay
!           cell(L)%deltaz = deltaz
!           ! Coordinates of the lower left corner of each cell
!           cell(L)%corner(1)  = -xmax + (I-1)*cell(L)%deltax 
!           cell(L)%corner(2)  = -ymax + (J-1)*cell(L)%deltay     
!           cell(L)%corner(3)  = -zmax + (K-1)*cell(L)%deltaz 
           ! Coordinates of the lower left corner of each cell
           cell(L)%corner(1)  = -xmax + (I-1)*cell_width
           cell(L)%corner(2)  = -ymax + (J-1)*cell_width     
           cell(L)%corner(3)  = -zmax + (K-1)*cell_width 
!           write(2,*) cell(L)%indexc, cell(L)%corner ! don't write if it is not necessary (computational very expensive)
           L = L + 1
        END DO
     END DO
  END DO

  

END SUBROUTINE setup_grid
