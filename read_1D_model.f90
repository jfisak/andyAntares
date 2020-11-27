SUBROUTINE read_1D_model() 

! Read 1D model data and allocet that data to the corresponding values of teh model grid cells

USE types

IMPLICIT NONE    

INTEGER                                   :: I, J, numbions, indexg, atom_number
! for reading from files
DOUBLE PRECISION                          :: junk
INTEGER                                   :: ios
INTEGER, PARAMETER                        :: maxrows = 6000000
DOUBLE PRECISION                          :: r, velo, dens, temp
DOUBLE PRECISION, DIMENSION(n_elements)   :: massfrac
DOUBLE PRECISION                          :: cell_index
CHARACTER(80)                             :: modelfile, jikrfile
! variables which are not needed in the code
!DOUBLE PRECISION                          :: delta_r, delta, delta2, tot_nd, tot_md
! (2) PoWR model
CHARACTER(100)                             :: powrfile
CHARACTER(100)                             :: line
DOUBLE PRECISION, PARAMETER                :: meanAtMass = 1.33
INTEGER                                         :: reading_grid
DOUBLE PRECISION, ALLOCATABLE                   :: boundaries(:)
DOUBLE PRECISION                                :: rPrev, rAct, width

modelfile=TRIM(inputmodelFile)

SELECT CASE (inputModel)
 CASE(0)
 OPEN (UNIT=11, FILE=modelfile)
 READ(11,*) T_eff
 READ(11,*) R_star
 ! READ(11,*) R_inf
 ! READ(11,*) V_inf
!  READ(11,*) M_dot
 ! READ(11,*) n_modelgrid
 
 
 add_mg = 1
  ! write(*,*) T_eff, R_star, R_inf, V_inf, M_dot, n_modelgrid
 
 ! R_star = R_star * r_sun
 ! WRITE(15, *) 0.D0, 0.D0, R_star/R_star
 ! write(*,*) '   ', T_eff, R_star, R_inf, V_inf, M_dot, n_modelgrid
 
 ! Allocate array for model grid structure.
 ! Cell n_modelgrid+1 is associated to propagation grid cells 
 ! which have no counterpart on the modelgrid  
  DO 
   READ(11, *, iostat = reading_grid) junk
   IF(reading_grid /= 0) EXIT
   n_modelgrid = n_modelgrid + 1
  END DO
  ALLOCATE (model_grid(n_modelgrid + add_mg))
  REWIND(11)
  READ(11,*) junk
  READ(11,*) junk
  DO I = 1, n_modelgrid
     ! Maybe better to calculate at the midle of the grid cell rather then at the outer boundary 
     READ(11,*) cell_index, r, velo, dens, temp!, massfrac
     model_grid(I)%rwind = r  * R_star
     model_grid(I)%vel = velo
     model_grid(I)%rho = dens
     model_grid(I)%T = temp ! should be temp 
     model_grid(I)%J = 0.D0 
     model_grid(I)%assoc_cells = 0
     ! write(*,'(A23, d14.5)') 'read_1D_model: rwind = ', r * R_star

     ALLOCATE (model_grid(I)%grid_comp(n_elements))
     DO J = 1, n_elements      
        numbions = elements(J)%nions
        ! write(*,*) 'read_1D_model: numbions = ', numbions
        ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
        atom_number = elements(J)%atom_number
        !model_grid(I)%grid_comp(J)%abund = massfrac(atom_number)        
        model_grid(I)%grid_comp(J)%abund = elements(J)%abundance
        !Calculate total number density for included species
        !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
        !model_grid(I)%grid_comp(J)%numb_den = tot_nd
     END DO
  END DO

  R_inf  = model_grid(n_modelgrid)%rwind
  V_inf  = model_grid(n_modelgrid)%vel
  ! temporary change
  V_inf = 30000.D+5
  write(*,*) 'read_1D_model: R_inf = ', R_inf/R_star, 'V_inf = ', V_inf

  ! setting properties


  ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
  ! All cells out of model grid set to 0 and associate to n_modelgrid. 
  ! Other cells will obtainde particular values with memory
  model_grid(n_modelgrid+1)%rwind = 0.D0
  model_grid(n_modelgrid+1)%vel   = 0.D0
  model_grid(n_modelgrid+1)%rho   = 0.D0     

 ! the model cell widths
 ! now it is a pont in the center of two neighbouring model cells
 ! the boundaries for the widths calculation
 ALLOCATE(boundaries(n_modelgrid + 1))
 boundaries(1) = R_star
 DO I = 1, n_modelgrid - 1
  boundaries(I + 1) = (model_grid(I)%rwind + model_grid(I + 1)%rwind)/2.0
 END DO
 boundaries(n_modelgrid + 1) = R_inf
 DO I = 1, n_modelgrid
  rPrev = boundaries(I)
  rAct = boundaries(I + 1)
  width = rAct - rPrev
  model_grid(I)%width = width
  ! write(*,*) 'read_1D_model: I = ', I, ' width = ', width
 END DO
 ! STOP 'read_1D_model: testing calculation of width'
 !______________________________________________________________________________________________
 ! (1) JIKR model
 !______________________________________________________________________________________________
 ! in this case we read input model from Jiri Krticka program...
 ! these files are in this form
 ! 1. number of row
 ! 2. wind radius
 ! 3. velocity
 ! 4. mass density
 ! 5. temperature
 ! 6. q (?)
 ! 7. mass loss rate
 CASE(1)
  write(99,*) 'we will read a model from Jiri Krticka program...'
  add_mg = 1
  CALL GET_ENVIRONMENT_VARIABLE("JIKRMODEL", jikrfile)
  IF(TRIM(jikrfile) == "") STOP "no input model file selected, &
                                   please set the variable JIKRMODEL"
  OPEN(UNIT=12,status='old',FILE=jikrfile)
   READ(12,*) T_eff, R_star, modelfile
  CLOSE(12)
  OPEN(UNIT=11,status='old',FILE=modelfile)
   ! at first the number of rows calculation...
   n_modelgrid = 0
  R_star = R_star * r_sun
!  T_eff = 37500
  DO I=1,maxrows
    READ(11,*,IOSTAT=ios) junk, junk, junk, junk, junk, junk, junk
   IF (ios /= 0) EXIT
   n_modelgrid = n_modelgrid + 1
  END DO
   ALLOCATE (model_grid(n_modelgrid + add_mg))
   REWIND(11)
  DO I=1,n_modelgrid
   READ(11,*) indexg, r, velo, dens, temp, junk, junk
     model_grid(I)%rwind = r  !  * R_star
     model_grid(I)%vel = velo * 1.E2
     model_grid(I)%rho = dens
     model_grid(I)%T = temp ! should be temp 
     model_grid(I)%J = 0.D0 
     model_grid(I)%assoc_cells = 0
     model_grid(I)%T = model_grid(I)%T / temp_factor
!     write(*,*) 'testing model grid...'
!     write(*,*) model_grid(I)%rwind, model_grid(I)%vel, &
!        model_grid(I)%rho, model_grid(I)%T, model_grid(I)%J, &
!        model_grid(I)%assoc_cells
     !IF (I /= 1) write(*,*) 'delta r: ', model_grid(I)%rwind - model_grid(I-1)%rwind
     ALLOCATE (model_grid(I)%grid_comp(n_elements))
     DO J = 1, n_elements      
        numbions = elements(J)%nions
        ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
        atom_number = elements(J)%atom_number
        model_grid(I)%grid_comp(J)%abund = elements(J)%abundance
        !Calculate total number density for included species
        !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
        !model_grid(I)%grid_comp(J)%numb_den = tot_nd
     END DO
   END DO
  CLOSE(11)
  IF(temp_factor /= 1.0) write(99,*) 'Warning, temperature structure is divided &
   by a temperature factor = ', temp_factor
  R_star = model_grid(1)%rwind
  R_inf  = model_grid(n_modelgrid)%rwind
  V_inf  = model_grid(n_modelgrid)%vel * 10.0**5
  ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
  ! All cells out of model grid set to 0 and associate to n_modelgrid. 
  ! Other cells will obtainde particular values with memory
  model_grid(n_modelgrid + 1)%rwind = 0.D0
  model_grid(n_modelgrid + 1)%vel   = 0.D0
  model_grid(n_modelgrid + 1)%rho   = 0.D0     
  ! calculating virtual particles from the selected input model
  !CALL virtual_particles(1)
 !______________________________________________________________________________________________
 ! (2) PoWR model
 !______________________________________________________________________________________________
 !
 ! 1. radius / R_*
 ! 2. radial velocity / km * s^(-1)
 ! 3. total numerical mass density
 ! 4. temperature / K
 CASE(2)
  write(99,*) 'we will read the TESTCASE from the PoWR code...'
  R_star = 20.066 * R_sun
  T_eff = 37000
  add_mg = 1
  CALL GET_ENVIRONMENT_VARIABLE("POWRMODEL", powrfile)
  IF(TRIM(powrfile) == "") STOP "no input model file selected, &
                                   please set the variable POWRMODEL"
  write(99,*) 'read_1D_model: powrfile = ', powrfile
  n_modelgrid = 0
  OPEN(UNIT=11, STATUS="old", FILE=TRIM(powrfile))
   DO
    read(11, *, IOSTAT=ios) line
    IF (ios /= 0) EXIT
    n_modelgrid = n_modelgrid + 1
   END DO
   ALLOCATE (model_grid(n_modelgrid + add_mg))
   write(99,*) 'read_1D_model: n_modelgrid = ', n_modelgrid
   REWIND(11)
   DO I=1,n_modelgrid
    READ(11,*) r, velo, dens, temp
    ! write(*,*) 'read_1D_model: I = ', I, ' r = ', r, ' velo = ', velo, ' dens = ', dens, ' temp = ', temp
    model_grid(I)%rwind = r * R_star
    model_grid(I)%vel = velo * 1.E5
    model_grid(I)%rho = dens * meanAtMass * mp_g
    model_grid(I)%T = temp 
    model_grid(I)%J = 0.D0 
    model_grid(I)%assoc_cells = 0
    ! print*, 'read_1D_model: testing model grid...'
    ! print*, 'read_1D_model: ', I, model_grid(I)%rwind, model_grid(I)%vel, &
    !  model_grid(I)%rho, model_grid(I)%T, model_grid(I)%J, &
    !  model_grid(I)%assoc_cells
    !IF (I /= 1) print*, 'delta r: ', model_grid(I)%rwind - model_grid(I-1)%rwind
    ALLOCATE (model_grid(I)%grid_comp(n_elements))
    DO J = 1, n_elements      
     numbions = elements(J)%nions
     ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
     atom_number = elements(J)%atom_number
     model_grid(I)%grid_comp(J)%abund = elements(J)%abundance
     !Calculate total number density for included species
     !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
     !model_grid(I)%grid_comp(J)%numb_den = tot_nd
    END DO
   END DO
  R_inf  = model_grid(1)%rwind
  write(*,*) 'read_1D_model: R_inf = ', R_inf / R_star
  V_inf  = model_grid(1)%vel
  write(99,*) 'read_1D_model: R_star = ', R_star, ' R_inf = ', R_inf
  ! Dummy cell to associate to propagation grid cells which have no
  ! representation on the model grid. All cells out of model grid
  ! set to 0 and associate to n_modelgrid. Other cells will obtainde
  ! particular values with memory
  model_grid(n_modelgrid + 1)%rwind = 0.D0
  model_grid(n_modelgrid + 1)%vel   = 0.D0
  model_grid(n_modelgrid + 1)%rho   = 0.D0     
  CLOSE(11)
  ! STOP 'read_1D_model: testing...'
 CASE DEFAULT
  write(99,*) 'the choice of the variable inputModel = ', inputModel, 'is not known...'
  STOP 'ENDING PROGRAM NOW...'
 END SELECT
  END SUBROUTINE read_1D_model
