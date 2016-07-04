SUBROUTINE main

  ! Propagate a bunch of photon packets through a stellar wind

! To run MPI (Message Passing Interface) it is necessary to include file mpif.h (Implicit Fortran MPI interfaces)
! to have output in an elegant manner, i.e. good formating
!#IFDEF MPI_ON
!  INCLUDE 'mpif.h' 
!#ENDIF

! Use module types (modul.f90)
  USE types


  IMPLICIT NONE

  INTEGER                           :: n_pack, iteration                                      !nx_cell, ny_cell, nz_cell
  INTEGER                           :: I, J, K, L, iseed, idx
  INTEGER, DIMENSION (9)            :: TT
  DOUBLE PRECISION, ALLOCATABLE     :: current_temp(:)
  CHARACTER(1)                      :: junk
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=2.D0/5.D0, lower_opa=0.01D0         ! Opacity for photons sent from the photosphere R_star = 10
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=2.D0/9.D0, lower_opa=0.1D0/18.D0    ! Opacity for photons sent from the photosphere R_star = 2
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=0.2, lower_opa=1.d0/200.d0          ! Opacity for photons sent from point sours
!  DOUBLE PRECISION                  :: xmax, ymax, zmax, deltax, deltay, deltaz 
!  DOUBLE PRECISION                  :: delta_cellx, delta_celly, delta_cellz 
!  DOUBLE PRECISION                  :: delta_opa, opa_cell

! Link data to identify program version
  CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
  COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
!  COMMON / RAN_SEED / idum


!  OPEN (UNIT=2, FILE='cells.dat')  
!  OPEN (UNIT=3, FILE='modelgrid.dat')
!  OPEN (UNIT=4, FILE='density.dat')

!  test = 0

!! Initialize MPI Parallelisation 
!#IFDEF MPI_ON
!  INTEGER  ierr  
!  ! Initialize the MPI execution environment 
!  CALL MPI_INIT(ierr)
!  ! If MPI routine completed successfully
!  IF (ierr .NE. MPI_SUCCESS) THEN
!    PRINT*,'Error starting MPI program. Terminating.'
!    ! If MPI routine failed
!    CALL MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
!  END IF
!  ! To determine rank within the set of processes
!  CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
!  ! To determine the number of processes
!  CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, ierr)
!  PRINT*, 'Number of tasks=',numtasks,' My rank=', my_rank
!#ENDIF


  ! Write Link Data (Program Version) to CPR file
  WRITE(*,'(2A)') '>>> Program started: Program Version from ',  &
                     LINK_DATE  !tag which returns the date the link was created as text
  WRITE(*,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)), &
           ' at host ', LINK_HOST(:IDX(LINK_HOST))

  ! Read input  
  CALL read_input(n_pack, iseed)
  PRINT*, 'read input'
  read*, junk

  ! Read composition
  CALL read_composition()

! Read atomic data (level information)
!  CALL read_atomic_data()
!   CALL read_atomic_characteristics()
  PRINT*, 'stop'

! Read transition data
  !CALL read_transitions()

!  PRINT*, 'stop'

!  STOP

  ! Initialing seed from the system time
  ! If we set iseed < 0 in input.dat then iseed will be randomly initializing from the system time
  ! otherwise, iseed will take a fix value given in the input file
  CALL DATE_AND_TIME(VALUES = TT)
  IF (iseed .LE. 0) THEN  
    iseed = TT(1)+70*(TT(2)+12*(TT(3)+31*(TT(5)+23*(TT(6)+59*TT(7)))))
  END IF

#ifdef MPI_ON
  ! For MPI parallel calculations each task needs its own random number seed
  ! This is achieved by adding an offset (in this case it is number 17, but 
  ! can be any other number) to the basic random number seed which depends 
  ! on the task's ID number. This has to be done for both "random" and 
  ! pre-defined seeds.
  iseed = iseed + my_rank*17
#endif


  ! The initial value of iseed (idum) should be set to different
  ! NEGATIVE integer values in order to obtain different random
  ! sequences. Seed is updated by ran2 once for each random number
  ! generated.
  idum = -iseed

  ! Only for debuging; if set the values > 0 then variou print out statement 
  ! will give information on a packet's history (depending on the actual value of debug)
  debug = 0

  ! Allocate array for photon packages.
  ! n_pack refers to a dummy package which can be used to sample
  ! packets properties while moved around.
  dummypackage = n_pack + 1
  ALLOCATE (package(n_pack + 1))

  ! Set up outflow (model grid)
  CALL setup_model_grid()
  xmax = xmax * R_star
  ymax = ymax * R_star 
  zmax = zmax * R_star
  ALLOCATE(current_temp(n_modelgrid+1))
  print*, 'model grid is set up'

  print*,'CHECK GRID SIZES'
  print*, xmax, R_inf, cell_width
  print*, xmax/R_star, R_inf/R_star, cell_width/R_star
  print*, (-xmax + nx_cell*cell_width)/R_star
  !STOP

  ! Set up of the propagation grid
  CALL setup_grid()
  CALL create_dynamical_grid_cell()
  print*, 'propagation grid is set up'

  L = 1
  DO I=1, nx_cell
     DO J=1, ny_cell
        DO K=1, nz_cell
           IF (I .EQ. nx_cell/2) THEN ! Only for one slice in the midle
              ! WRITE(3, *) cell(L)%corner, model_grid(cell(L)%model_index)%rho
              ! WRITE(4, *) model_grid(cell(L)%model_index)%rwind, model_grid(cell(L)%model_index)%rho
              ! WRITE(4, *) model_grid(cell(L)%model_index)%rho
           END IF
           L = L + 1
        END DO
     END DO
  END DO         
  print*, 'Check model grid done'

  ! Checking if the analitic solution for the escape probability
  ! (e^(-tau)) is in agreement with the calculated one using our
  ! propagation procedure

  ! Increament of opacity (for testing)
  !  delta_opa = (upper_opa - lower_opa)/nopa

  ! Update model grid properties (model will be updated after 
  ! consistance temperature calculation from teh radiation field)

  current_temp = 0.D0

  DO iteration = 1, 1

     CALL update_grid(iteration)
     PRINT*, 'Update grid finished' 
!     print*, 'model_grid(:)%T = ', model_grid(:)%T
     IF ( MAXVAL((model_grid(:)%T - current_temp) / model_grid(:)%T) .LT. 0.05D0) EXIT

     current_temp = model_grid(:)%T
!     PRINT*, 'iteration:', iteration, current_temp

     ! Loop over the numer of different opacity (nopa)
     ! DO I=1,nopa
     ! print*, 'tau loop',  I
     ! Opacity for tau calculation
     ! opa_cell = lower_opa  +  (I-1)*delta_opa
     ! print*,   opa_cell * (xmax-R_star)
     ! print*, R_star
     ! stop

     ! Initialisation of photon packages from the photosphere
     CALL init_photsphere(n_pack) 
     print*, 'photons initialised'

     ! Initalisation of photon packages from point sourse
     ! CALL init_photonpack(n_pack)

     ! Propagation of the photon in 3D grid
     ! CALL propagation(n_pack, opa_cell, lower_opa, delta_opa)
     print*, 'update packages'
     CALL update_packages(n_pack)
     print*, 'Number of destoyed packages =', destroyed_pack

  END DO 

  IF (iteration .GE. 10) print*, 'No convergency'
 
  print*, 'do spectrum'
  CALL do_spectrum(n_pack)
  print*, 'do finalize'

     ! END DO
    

!#ifdef MPI_ON
!  call MPI_FINALIZE(ierr)
!#endif

END SUBROUTINE main
