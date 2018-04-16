SUBROUTINE main

  ! Propagate a bunch of photon packets through a stellar wind

! Use module types (modul.f90)
  USE types


  IMPLICIT NONE

  INTEGER                           :: n_pack, iteration!, nx_cell, ny_cell, nz_cell
  INTEGER                           :: I, J, K, L, iseed, idx
  INTEGER, DIMENSION (9)            :: TT
  DOUBLE PRECISION, ALLOCATABLE     :: current_temp(:)
  CHARACTER(1)                      :: junk
  CHARACTER                         :: n_dummy_packs_char
  CHARACTER(2)                      :: chnum_threads
  INTEGER                           :: num_threads, stat
! parallelized part
! definition of MPI variables
! #if mpi==1
 INTEGER                            :: numtasks
 INTEGER                            :: rc
 CHARACTER(20)                      :: outputfile
! #endif
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=2.D0/5.D0, lower_opa=0.01D0         ! Opacity for photons sent from the photosphere R_star = 10
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=2.D0/9.D0, lower_opa=0.1D0/18.D0    ! Opacity for photons sent from the photosphere R_star = 2
!  DOUBLE PRECISION, PARAMETER       :: upper_opa=0.2, lower_opa=1.d0/200.d0          ! Opacity for photons sent from point sours
!  DOUBLE PRECISION                  :: xmax, ymax, zmax, deltax, deltay, deltaz 
!  DOUBLE PRECISION                  :: delta_cellx, delta_celly, delta_cellz 
!  DOUBLE PRECISION                  :: delta_opa, opa_cell
INTEGER                              :: nphit, loc_n_pack
INTEGER, PARAMETER                   :: max_packs = 1e7

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
 ! write(*,*) 'mpi initialization', ' my_rank = ', my_rank, &
 !  ' n_tasks = ', n_tasks
 ! write(*,*) 'main: my_rank = ', my_rank
 ! write(*,*) outputfile
 ! write(chnum_threads, "(A)"), my_rank
 ! outputfile = '/packets.'//chnum_threads
 ! write(*,*) 'main: outputfile = ', outputfile
 ! write(*,*) 'main: outputfile = ', outputfile
#endif
write(outputfile,"(A6,I3.3)") "output", my_rank
OPEN(99, FILE=outputfile) 
write(*,*) 'mpi initialization: my_rank = ', my_rank, &
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
 ! Read composition
 write(99,*) 'read_composition'
 CALL read_composition()

 ! Initialing seed from the system time
 ! If we set iseed < 0 in input.dat then iseed
 ! will be randomly initializing from the system time
 ! otherwise, iseed will take a fix value given in
 ! the input file
 CALL DATE_AND_TIME(VALUES = TT)
 IF (iseed .LE. 0) THEN  
  iseed = TT(1)+70*(TT(2)+12*(TT(3)+31*(TT(5)+23*(TT(6)+59*TT(7)))))
 END IF


 ! The initial value of iseed (idum) should be set to different
 ! NEGATIVE integer values in order to obtain different random
 ! sequences. Seed is updated by ran2 once for each random number
 ! generated.
 idum = -iseed

! Only for debuging; if set the values > 0 then variou print out statement 
! will give information on a packet's history (depending on the actual value of debug)
debug = 0

! Allocate array for photon packages.
ALLOCATE (package(n_pack + 1))

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
 
ALLOCATE(current_temp(n_modelgrid + add_mg))
write(99,*) 'model grid is set up'
write(99,*) 'setup propagation grid'

! Set up of the propagation grid
CALL setup_grid2()
! connects the propagation grid with the model grid
write(99,*) 'propagation grid is set up'
CALL connection_prop_model_grid()


current_temp = 0.D0
iteration = 0

! iteration = 0
! OPEN(20, FILE='temp_structure.dat')
DO iteration = 1,1
 ! definition of counters
 ! iteration = iteration + 1
 IF (iteration .GE. 10) print*, 'No convergency'
 CALL update_grid(iteration)
 IF(iteration == 100) STOP 'too many iteration in the subroutine main'
 write(*,*) 'Update grid finished' 
! WRITE(20,*) '# ITERATION: ', iteration
 DO I = 1, n_modelgrid
! WRITE(20,*) I, model_grid(I)%T
 END DO
!   print*, 'model_grid(:)%T = ', model_grid(:)%T
!   IF ( MAXVAL(abs(model_grid(:)%T - current_temp) / model_grid(:)%T) .LT. 0.05D-1) EXIT
  current_temp = model_grid(:)%T
!  PRINT*, 'iteration:', iteration, current_temp

! Loop over the numer of different opacity (nopa)
! DO I=1,nopa
! print*, 'tau loop',  I
! Opacity for tau calculation
! opa_cell = lower_opa  +  (I-1)*delta_opa
! print*,   opa_cell * (xmax-R_star)
! print*, R_star
! stop
! n_pack will compute in several loops to save some memory
 nphit = INT(n_pack / max_packs)
! do 03
 nphit = 1
!  DO J = 1, nphit
! Initialisation of photon packages from the photosphere
 CALL init_photsphere(n_pack) 
 write(*,*) 'photons initialised'

! Initalisation of photon packages from point source
! CALL init_photonpack(n_pack)

! Propagation of the photon in 3D grid
! CALL propagation(n_pack, opa_cell, lower_opa, delta_opa)
 write(*,*) 'update packages'
 CALL update_packages(n_pack)
 write(*,*) 'Number of destoyed packages =', destroyed_pack
! end do 03
!  END DO
END DO 
! CLOSE(20)

 
 print*, 'do spectrum'
 ! CALL do_spectrum(n_pack)
 print*, 'do finalize'
 ! it will save some important output
 CALL save_output(1)
 CALL save_output(2)
 ! temp structure and occupation numbers
! CALL save_output(3)
 CALL save_output(4)

#ifdef MPI_ON
 call MPI_FINALIZE(ierr)
#endif
CLOSE(2)

END SUBROUTINE main
