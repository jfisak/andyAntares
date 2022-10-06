SUBROUTINE main

  ! Propagate a bunch of photon packets through a stellar wind

! Use module types (modul.f90)
  USE types


  IMPLICIT NONE

  INTEGER                           :: n_pack, iteration!, nx_cell, ny_cell, nz_cell
  INTEGER                           :: iseed, idx
  INTEGER, DIMENSION (9)            :: TT
  DOUBLE PRECISION, ALLOCATABLE     :: current_temp(:)
  INTEGER                               :: cur_parameter=0
  REAL                                  :: time0_agcreation, time1_agcreation
  REAL                                  :: time0_agconnwpg, time1_agconnwpg
  REAL                                  :: time0_pp, time1_pp
  REAL                                  :: time_cre, time_con, time_pp
! parallelized part
! definition of MPI variables
INTEGER                              :: nphit

LOGICAL                                 :: propmod_file_exists
CHARACTER(60)                           :: propmod_file

LOGICAL                                 :: timing = .true.


! Link data to identify program version
CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
! COMMON / RAN_SEED / idum

! the Saha constant calculation
saha_const = 5.D-1 * (h**2/(2.0*pi*me_g*BOLK))**1.5
! CALL EXECUTE_COMMAND_LINE('figlet "3D WIND CODE"')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MPI INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
my_rank = 0
#if mpi==1
 CALL MPI_INIT(ierr)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, ierr)
#endif


 CALL save_output(cur_parameter)
 write(99,*) 'mpi initialization: my_rank = ', my_rank, &
  ' n_tasks = ', n_tasks
 WRITE(99,'(2A)') '>>> Program started: Program Version from ',  &
  LINK_DATE ! tag which returns the date the link was created as text
 WRITE(99,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)), &
  ' at host ', LINK_HOST(:IDX(LINK_HOST))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! READ INPUT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(99,*) 'read input'
 CALL read_input(n_pack, iseed)
 CALL analyse_input()
 ! Allocate array for photon packages.
 ALLOCATE (package(n_pack + 2))
 ! write(*,*) 'main: |package| = ', SIZE(package)

 CALL find_unfinished_run()
 ! Read composition
 write(99,*) 'read_composition'
 CALL read_composition()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MONTE CARLO SEED INITIALIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Initialing seed from the system time
 ! If we set iseed < 0 in input.dat then iseed
 ! will be randomly initializing from the system time
 ! otherwise, iseed will take a fix value given in
 ! the input file
 IF(iseed > 0 .AND. n_tasks > 1) STOP 'iseed > 0 & n_tasks > 1'
 CALL DATE_AND_TIME(VALUES = TT)
 IF (iseed .LE. 0) THEN  
  iseed = TT(1)+70*(TT(2)+12*(TT(3)+31*(TT(5)+23*(TT(6)+59*TT(7)))))
 END IF
 iseed = iseed + my_rank * 17 ! add + threadID * primeNumber


 ! The initial value of iseed (idum) should be set to different
 ! NEGATIVE integer values in order to obtain different random
 ! sequences. Seed is updated by ran2 once for each random number
 ! generated.
 idum = -iseed
 write(99,*) 'main: iseed = ', iseed, ' idum = ', idum

! Only for debuging; if set the values > 0 then variou print out statement 
! will give information on a packet's history (depending on the actual value of debug)
! 0 -- no debug mode
! 2 -- packet propagation debugging
debug = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETTING OF PROPAGATION AND MODEL GRID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(propmod_file,"(A, A12)") TRIM(outputfolder), '/propmod.dat'
INQUIRE(FILE=propmod_file, EXIST=propmod_file_exists)

IF(saved_grid == 1 .and. propmod_file_exists) THEN
 CALL read_propmod_grid()
ELSE
 ! write(*,*) 'setting up model grid'
 CALL setup_model_grid()
 
 ! if model_type == 3 xyzmax are already calculated in setup_model_grid
 IF(model_type /= 3) THEN
   xmax = R_inf + 0.5 * R_sun
   ymax = R_inf + 0.5 * R_sun 
  IF(model_type == 1) THEN
   zmax = R_inf + 0.5 * R_sun
  ELSE IF (model_type == 2 .AND. inputmodel == 1) THEN
   zmax = Z_inf! + R_sun
  ELSE IF (model_type == 2 .AND. inputmodel == 2) THEN
   zmax = R_inf + 0.5 * R_sun
  ELSE
   STOP 'main: non-known model type'
  END IF
 END IF
 
 write(99,*) 'model grid is set up'
 write(99,*) 'setup propagation grid'
 write(99,*) 'xmax = ', xmax/R_star, ' ymax = ', ymax/R_star, ' zmax = ', zmax/R_star
  
 
 ! Set up of the propagation grid
 write(99,*) 'setting up the propagation grid'
 ! write(*,*) 'setting up the propagation grid'
 CALL setup_propgrid()
 
 ! write(*,*) 'connecting prop and mod grids'
 CALL connection_prop_model_grid()
 ! write(*,*) 'prop and mod grids are connected'
END IF ! saved propmod grid

! save propmod_grid?
IF(saved_grid == 1 .and. .not. propmod_file_exists) THEN
#if mpi==1
 IF(my_rank == 0) THEN
#endif
 CALL save_propmod_grid()
#if mpi==1
 END IF
#endif
END IF
! connects the propagation grid with the model grid
write(99,*) 'propagation grid is set up'

 ! IF(my_rank == 0) THEN
 !  CALL save_output(8)
 !  ! STOP 'testing the propagation grid'
 ! END IF





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THE MAIN ITERATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(current_temp(n_modelgrid + add_mg))
current_temp = 0.D0
iteration = 0

! OPEN(20, FILE='temp_structure.dat')
DO iteration = 1,1
 ! definition of counters
 ! iteration = iteration + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETTING THE PLASMA STATE IN THE MODEL GRID!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 IF (iteration .GE. 10) write(99,*) 'No convergency'
 ! write(*,*) 'updating modGrid'
 CALL update_grid(iteration)
 ! write(*,*) 'modGrid was updated'
 ! if (iteration == 1 .AND. inputpopfile .NE. '') then
 !  CALL read_populations()
 ! end if
 CALL i_ion_recomb(1)
 IF(iteration == 100) STOP 'too many iteration in the subroutine main'
 write(99,*) 'Update grid finished' 
! WRITE(20,*) '# ITERATION: ', iteration
  current_temp = model_grid(:)%T
! do 03
 nphit = 1
!  DO J = 1, nphit
! Initialisation of photon packages from the photosphere
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PACKET MACHINERY !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! write(*,*) 'initialization of photosphere'
 CALL init_photsphere(n_pack) 
 ! write(*,*) 'photosphere was initialized'

! Propagation of the photon in 3D grid
! CALL propagation(n_pack, opa_cell, lower_opa, delta_opa)
 write(99,*) 'update packages'
 ! write(*,*) 'update packages'
 IF(timing) THEN
  CALL cpu_time(time0_pp)
 END IF
 CALL update_packages(n_pack)
 IF(timing) THEN
  CALL cpu_time(time1_pp)
  write(99,*) 'time0_pp = ', time0_pp, ' time1_pp = ', time1_pp
 END IF
 ! write(*,*) 'packets were updated'
 write(99,*) 'Number of destoyed packages =', destroyed_pack
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DISTRIBUTE VARIABLES AMONG DIFFERENT TASKS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#if mpi==1
 CALL mpi_distribute_estimators()
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
END DO ! iteration (now of temperature structure)
! CLOSE(20)

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SAVE OUTPUT FILES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 write(99,*) 'do finalize'
 ! it will save some important output
 CALL save_output(1)
 CALL save_output(2)
 ! temp structure and occupation numbers
 IF(my_rank == 0) THEN
  CALL save_output(3)
  CALL save_output(7)
  ! partition function
  CALL save_output(9)
 END IF
 CALL save_output(5)
 ! info o gridu
 ! IF(my_rank == 0) CALL save_output(6)
 ! erase all temporary files with photons

CLOSE(2)
CLOSE(99)

END SUBROUTINE main
