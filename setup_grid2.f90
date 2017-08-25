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
  INTEGER                                :: xp, xm, yp, ym, zp, zm
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
  basic_cell_width(1) = 2.E0 * xmax / DBLE(nx_cell)
  basic_cell_width(2) = 2.E0 * ymax / DBLE(ny_cell)
  basic_cell_width(3) = 2.E0 * zmax / DBLE(nz_cell)

  L = 1
  DO I=1, nx_cell
   DO J=1, ny_cell
    DO K=1, nz_cell
     ! Index(number) of each cell in x,y, and z direction
     !dyn_cell(L)%indexc(1) = I
     !dyn_cell(L)%indexc(2) = J
     !dyn_cell(L)%indexc(3) = K
     ! Coordinates of the lower left corner of each cell
     dyn_cell(L)%corner(1)  = - xmax + DBLE((I - 1)) * basic_cell_width(1)
     dyn_cell(L)%corner(2)  = - ymax + DBLE((J - 1)) * basic_cell_width(2)     
     dyn_cell(L)%corner(3)  = - zmax + DBLE((K - 1)) * basic_cell_width(3) 
     ! cell width
     dyn_cell(L)%width(1) = basic_cell_width(1)
     dyn_cell(L)%width(2) = basic_cell_width(2)
     dyn_cell(L)%width(3) = basic_cell_width(3)
     ! number of down cell is equal to zero
     dyn_cell(L)%down_cell = 0
     dyn_cell(L)%up_cell = 0
     ! calculation of neighbors
     ! x+
     IF(I == nx_cell) THEN
      xp = -99
     ELSE
      xp = L + ny_cell * nz_cell
     END IF
     ! x-
     IF(I == 1) THEN
      xm = -99
     ELSE
      xm = L - ny_cell * nz_cell
     END IF
     ! y+
     IF(J == ny_cell) THEN
      yp = -99
     ELSE
      yp = L + nz_cell
     END IF
     ! y-
     IF(J == 1) THEN
      ym = -99
     ELSE
      ym = L - nz_cell
     END IF
     ! z+
     IF(K == nz_cell) THEN
      zp = -99
     ELSE
      zp = L + 1
     END IF
     ! z-
     IF(K == 1) THEN
      zm = -99
     ELSE
      zm = L - 1
     END IF
     ! 
     dyn_cell(L)%neighbor(1) = xp
     dyn_cell(L)%neighbor(2) = xm
     dyn_cell(L)%neighbor(3) = yp
     dyn_cell(L)%neighbor(4) = ym
     dyn_cell(L)%neighbor(5) = zp
     dyn_cell(L)%neighbor(6) = zm
     L = L + 1
    END DO
   END DO
  END DO
  ! the maximal number of cells is now equal to L
  max_n_dcell = Ngrid

 IF(dyngrid /= 0) THEN
  DO I = 1, Ngrid
   CALL create_dynamical_grid_cells(I, max_n_dcell)
   !print*, 'max_n_dcell = ', max_n_dcell
  END DO
 END IF
  !OPEN(15,FILE='dyn_cells.dat')
  ! DO I = 1, max_n_dcell
  !  write(15,*) dyn_cell(I)%corner, dyn_cell(I)%width, dyn_cell(I)%neighbor
  ! END DO
  !CLOSE(15)
! DO I = 1, max_n_dcell
!  IF(dyn_cell(I)%up_cell == 0) THEN
!   ! what is the basic cell of this dynamical cell
!   actCell = I
!   DO WHILE(dyn_cell(actCell)%down_cell /= 0)
!    actCell = dyn_cell(actCell)%down_cell
!   END DO
!   basicCell = actCell
!  END IF
!
! END DO

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
