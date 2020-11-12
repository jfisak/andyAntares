SUBROUTINE main

  ! Propagate a bunch of photon packets through a stellar wind

! Use module types (modul.f90)
  USE types


  IMPLICIT NONE

  INTEGER                           :: n_pack, iteration!, nx_cell, ny_cell, nz_cell
  INTEGER                           :: iseed, idx
  INTEGER, DIMENSION (9)            :: TT
  DOUBLE PRECISION, ALLOCATABLE     :: current_temp(:)
! parallelized part
! definition of MPI variables
INTEGER                              :: nphit

! Link data to identify program version
CHARACTER LINK_DATE*30, LINK_USER*10, LINK_HOST*60
COMMON / COM_LINKINFO / LINK_DATE, LINK_USER, LINK_HOST
! COMMON / RAN_SEED / idum

! CALL EXECUTE_COMMAND_LINE('figlet "3D WIND CODE"')
my_rank = 0
#if mpi==1
 CALL MPI_INIT(ierr)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, n_tasks, ierr)
#endif


 CALL save_output(0)
 write(99,*) 'mpi initialization: my_rank = ', my_rank, &
  ' n_tasks = ', n_tasks
! CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
! STOP
 WRITE(99,'(2A)') '>>> Program started: Program Version from ',  &
  LINK_DATE ! tag which returns the date the link was created as text
 WRITE(99,'(4A)') '>>> created by ', LINK_USER(:IDX(LINK_USER)), &
  ' at host ', LINK_HOST(:IDX(LINK_HOST))

 ! Read input  
 write(99,*) 'read input'
 CALL read_input(n_pack, iseed)
 ! Allocate array for photon packages.
 ALLOCATE (package(n_pack + 2))
 ! write(*,*) 'main: |package| = ', SIZE(package)

 CALL find_unfinished_run()
 ! Read composition
 write(99,*) 'read_composition'
 CALL read_composition()

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
debug = 0

! Set up outflow (model grid)
CALL setup_model_grid()
! create virtual particles for the given model cell
IF (dyngrid /= 0) CALL virtual_particles(model_type)
 xmax = R_inf! + R_sun
 ymax = R_inf! + R_sun 
IF(model_type == 1) THEN
 zmax = R_inf! + R_sun
ELSE IF (model_type == 2 .AND. inputmodel == 1) THEN
 zmax = Z_inf! + R_sun
ELSE
 STOP 'main: non-known model type'
END IF
write(99,*) 'model grid is set up'
write(99,*) 'setup propagation grid'
write(99,*) 'xmax = ', xmax/R_star, ' ymax = ', ymax/R_star, ' zmax = ', zmax/R_star
 
ALLOCATE(current_temp(n_modelgrid + add_mg))

CALL DATE_AND_TIME(VALUES = TT)
write(*,*) TT(5), TT(6), TT(7), TT(8)
! Set up of the propagation grid
CALL setup_grid2()
CALL DATE_AND_TIME(VALUES = TT)
write(*,*) TT(5), TT(6), TT(7), TT(8)
! connects the propagation grid with the model grid
write(99,*) 'propagation grid is set up'
CALL connection_prop_model_grid()
CALL DATE_AND_TIME(VALUES = TT)
write(*,*) 'connected prop/model grid', TT(5), TT(6), TT(7), TT(8)

 IF(my_rank == 0) THEN
  CALL save_output(8)
  CALL DATE_AND_TIME(VALUES = TT)
  write(*,*) TT(5), TT(6), TT(7), TT(8)
  ! STOP 'testing the propagation grid'
 END IF
current_temp = 0.D0
iteration = 0

! OPEN(20, FILE='temp_structure.dat')
DO iteration = 1,1
 ! definition of counters
 ! iteration = iteration + 1
 IF (iteration .GE. 10) write(99,*) 'No convergency'
 CALL update_grid(iteration)
 CALL i_ion_recomb(1)
 IF(iteration == 100) STOP 'too many iteration in the subroutine main'
 write(99,*) 'Update grid finished' 
! WRITE(20,*) '# ITERATION: ', iteration
  current_temp = model_grid(:)%T
! do 03
 nphit = 1
!  DO J = 1, nphit
! Initialisation of photon packages from the photosphere
 CALL init_photsphere(n_pack) 

! Initalisation of photon packages from point source
! CALL init_photonpack(n_pack)

! Propagation of the photon in 3D grid
! CALL propagation(n_pack, opa_cell, lower_opa, delta_opa)
 write(99,*) 'update packages'
 CALL update_packages(n_pack)
 write(99,*) 'Number of destoyed packages =', destroyed_pack
 CALL DATE_AND_TIME(VALUES = TT)
 write(*,*) 'packet prop: ', TT(5), TT(6), TT(7), TT(8)
#if mpi==1
 CALL mpi_distribute_estimators()
 CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
END DO ! iteration (now of temperature structure)
! CLOSE(20)

 
 ! write(99,*) 'do spectrum'
 ! CALL do_spectrum(n_pack)
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
! #if mpi==1
!  ! will delete all temporary files
!  IF(my_rank == 0) THEN
!   cmdcommand = 'rm '//TRIM(outputfolder)//'/'//TRIM(temp_filename)//'*.dat'
!  ! write(*,*) 'main: cmdcommand = ', cmdcommand
!  CALL SYSTEM(cmdcommand)
!  END IF
!  CALL MPI_FINALIZE(ierr)
! #else
!  cmdcommand = 'rm '//TRIM(outputfolder)//'/'//TRIM(temp_filename)//'*.dat'
!  write(*,*) 'main: cmdcommand = ', cmdcommand
!  CALL SYSTEM(cmdcommand)
! #endif

END SUBROUTINE main
