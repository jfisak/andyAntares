  SUBROUTINE setup_modelgrid() 

  USE types

  IMPLICIT NONE    

  INTEGER                           :: I, J, K, L
  DOUBLE PRECISION, DIMENSION(3)    :: cell_center
  DOUBLE PRECISION                  :: M_dot, r

  M_dot = 5.D0

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
!          Absolute radius of the grid cell 
           r =SQRT( (cell(L)%corner(1) + cell_width/2.D0)**2 
              + (cell(L)%corner(2) + cell_width/2.D0))**2    &
              + (cell(L)%corner(3) + cell_width/2.D0))**2)
           L = L + 1
        END DO
     END DO
  END DO

  END SUBROUTINE setup_modelgrid
