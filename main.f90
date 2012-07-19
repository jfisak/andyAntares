SUBROUTINE main

  USE types

! Purpose: to calculet the escape probability 

  IMPLICIT NONE

  INTEGER                           :: n_pack, n_bin!, nx_cell, ny_cell, nz_cell
  INTEGER                           :: I, J, K, L, iseed, idx

  INTEGER, PARAMETER                :: Nmax = 1000000000, nopa =20
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=2.D0/5.D0, lower_opa=0.01D0    ! Opacity for photons sent from the photosphere R_star = 10
  DOUBLE PRECISION, PARAMETER       :: upper_opa=2.D0/9.D0, lower_opa=0.1D0/18.D0    ! Opacity for photons sent from the photosphere    R_star = 2
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=0.2, lower_opa=1.d0/200.d0         ! Opacity for photons sent from point sours
  INTEGER, DIMENSION (9)            :: TT
!  DOUBLE PRECISION                  :: xmax, ymax, zmax, deltax, deltay, deltaz, 
!  DOUBLE PRECISION                  :: delta_cellx, delta_celly, delta_cellz, 
  DOUBLE PRECISION                  :: delta_opa, opa_cell

! Link data to identify program version
  CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
  COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
!  COMMON / RAN_SEED / idum

  OPEN (UNIT=2, FILE='cells.dat')  
  OPEN (UNIT=3, FILE='modelgrid.dat')
!  OPEN (UNIT=4, FILE='density.dat')

!  test = 0

! Write Link Data (Program Version) to CPR file
  WRITE(*,'(2A)') '>>> Program started: Program Version from ',  &
                     LINK_DATE
  WRITE(*,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)), &
           ' at host ', LINK_HOST(:IDX(LINK_HOST))
  CALL read_input(n_pack, n_bin, iseed)

! Initialing seed from the system time
! If we set iseed < 0 in input.dat then iseed will be initializing from the system time
! otherwise, iseed will take a value given in the input file
  print*, 'read input'
  CALL DATE_AND_TIME(VALUES = TT)
  IF (iseed .LE. 0) THEN  
    iseed = TT(1)+70*(TT(2)+12*(TT(3)+31*(TT(5)+23*(TT(6)+59*TT(7)))))
  END IF

! The initial value of iseed (idum) should be set to different NEGATIVE integer values 
! in order to obtain different random sequences. Seed is updated by ran2 once for each 
! random number generated.
  idum = -iseed

  debug = 0

!  R_star = R_star * r_sun

! Number of grid cells
  Ngrid = nx_cell * ny_cell * nz_cell

  IF (Ngrid .GT. Nmax) THEN
     print*, 'ERROR: N > Nmax', Ngrid
     STOP 
  END IF

! Define length of dinamic arrays (cell and package)
  ALLOCATE (cell(Ngrid))
  ALLOCATE (package(n_pack))
  ALLOCATE (model_grid(n_modelgrid))

!  TYPE(grid_cell), DIMENSION(N) :: cell  

! Size of the grid cells in x,y, and z direction (now they are with the same size i.e. regular gred)
  cell_width = 2.D0*xmax/nx_cell
!  print*, cell_width
!! Size of the grid cells in x,y, and z direction (now they are with the same size i.e. regular gred)
!  deltax = 2.D0*xmax/nx_cell
!  deltay = 2.D0*ymax/ny_cellR
!  deltaz = 2.D0*zmax/nz_cell

! Set up of the gred
  CALL setup_grid()

!  L = 1
!  DO I=1, nx_cell
!     DO J=1, ny_cell
!        DO K=1, nz_cell
!           ! Index(number) of each cell in x,y, and z direction
!           cell(L)%indexc(1) = I
!           cell(L)%indexc(2) = J
!           cell(L)%indexc(3) = K
!           ! Size of each grid cells in x,y, and z direction
!           cell(L)%deltax = deltax
!           cell(L)%deltay = deltay
!           cell(L)%deltaz = deltaz
!           ! Coordinates of the lower left corner of each cell
!           cell(L)%corner(1)  = -xmax + (I-1)*cell(L)%deltax 
!           cell(L)%corner(2)  = -ymax + (J-1)*cell(L)%deltay     
!           cell(L)%corner(3)  = -zmax + (K-1)*cell(L)%deltaz 
!           ! Coordinates of the lower left corner of each cell
!           cell(L)%corner(1)  = -xmax + (I-1)*cell_width
!           cell(L)%corner(2)  = -ymax + (J-1)*cell_width     
!           cell(L)%corner(3)  = -zmax + (K-1)*cell_width 
!           write(2,*) cell(L)%indexc, cell(L)%corner ! don't write if it is not necessary (computational very expensive)
!           L = L + 1
!        END DO
!     END DO
!  END DO


  CALL toy_model()     ! Set up outflow

  L = 1
  DO I=1, nx_cell
     DO J=1, ny_cell
        DO K=1, nz_cell
           IF (I .EQ. nx_cell/2) THEN ! Only for one slice in the midle
               WRITE(3, *) cell(L)%corner, model_grid(cell(L)%model_index)%rho
!               WRITE(4, *) model_grid(cell(L)%model_index)%rho
           END IF
           L = L + 1
        END DO
     END DO
  END DO         


! Checking if the analitic solution for the escape probability (e^(-tau)) is in agreement with the
! calculated one using our propagation procedure

! Increament of opacity (for testing)
  delta_opa = (upper_opa - lower_opa)/nopa

! Loop over the numer of different opacity (nopa)
  DO I=1,nopa
     print*, 'tau loop',  I
!    Opacity for tau calculation
     opa_cell = lower_opa  +  (I-1)*delta_opa
!     print*,   opa_cell * (xmax-R_star)
!     print*, R_star
!     stop
!    Initalisation of photon packages from the photosphere
     CALL init_photsphere(n_pack) 
!    Initalisation of photon packages from point sourse
!     CALL init_photonpack(n_pack)
!    Propagation of the photon in 3D gred
     CALL propagation(n_pack, opa_cell, lower_opa, delta_opa)
  END DO
    
END SUBROUTINE main
