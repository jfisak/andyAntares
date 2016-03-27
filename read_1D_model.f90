  SUBROUTINE read_1D_model() 

  ! Read 1D model data and allocet that data to the corresponding values of teh model grid cells

  USE types

  IMPLICIT NONE    

  INTEGER                                   :: I, J, M, numbions, indexg, atom_number
  INTEGER                                   :: ios
  INTEGER, PARAMETER                        :: maxrows = 6000000
  DOUBLE PRECISION                          :: delta_r, delta, delta2, tot_nd, tot_md
  DOUBLE PRECISION                          :: r, velo, dens, temp, junk
  DOUBLE PRECISION, DIMENSION(n_elements)   :: massfrac
  CHARACTER(20)                             :: modelfile 

 SELECT CASE (inputModel)
  CASE(0)
  OPEN (UNIT=11, FILE='model_data.dat')

  READ(11,*) T_eff
  READ(11,*) R_star
!  READ(11,*) R_inf
!  READ(11,*) V_inf
!  READ(11,*) M_dot
  READ(11,*) n_modelgrid

!  print*, T_eff, R_star, R_inf, V_inf, M_dot, n_modelgrid

  R_star = R_star * r_sun
  ! WRITE(15, *) 0.D0, 0.D0, R_star/R_star

!  R_inf  = R_inf  * R_star
!  V_inf  = V_inf  * 1.D5
!  M_dot  = M_dot  * m_sun / (3600.D0*24.D0*365.25D0)

!  print*, '   ', T_eff, R_star, R_inf, V_inf, M_dot, n_modelgrid

 ! Allocate array for model grid structure.
  ! Cell n_modelgrid+1 is associated to propagation grid cells 
  ! which have no counterpart on the modelgrid  
  ALLOCATE (model_grid(n_modelgrid + 1))

  DO I = 1, n_modelgrid
     ! Maybe better to calculate at the midle of the grid cell rather then at the outer boundary 
     READ(11,*) indexg, r, velo, dens, temp, massfrac
     model_grid(I)%rwind = r  * R_star
     model_grid(I)%vel = velo * 1.D5
     model_grid(I)%rho = dens
     model_grid(I)%T = 8000. ! should be temp 
     model_grid(I)%J = 0.D0 
     model_grid(I)%assoc_cells = 0
     !Total mass density of grid cell I
!     tot_md = M_dot / (4.D0 * pi * (model_grid(I)%rwind)**2 * model_grid(I)%vel)     
     !WRITE(15, *) 0.D0, 0.D0, model_grid(I)%rwind/R_star, model_grid(I)%rho
!     print*, I, model_grid(I)%rwind, model_grid(I)%vel, model_grid(I)%rho, model_grid(I)%T
!     print*, I, model_grid(I)%rwind, model_grid(I)%vel, model_grid(I)%rho, tot_md

     ALLOCATE (model_grid(I)%grid_comp(n_elements))
     DO J = 1, n_elements      
        numbions = elements(J)%nions
        ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
        atom_number = elements(J)%atom_number
        model_grid(I)%grid_comp(J)%abund = massfrac(atom_number)        
        !Calculate total number density for included species
        !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
        !model_grid(I)%grid_comp(J)%numb_den = tot_nd
     END DO
  END DO

  R_inf  = model_grid(n_modelgrid)%rwind
  V_inf  = model_grid(n_modelgrid)%vel


  ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
  ! All cells out of model grid set to 0 and associate to n_modelgrid. 
  ! Other cells will obtainde particular values with memory
  model_grid(n_modelgrid+1)%rwind = 0.D0
  model_grid(n_modelgrid+1)%vel   = 0.D0
  model_grid(n_modelgrid+1)%rho   = 0.D0     


! Only for testing 
! DO I=1, Ngrid
!   r =SQRT( (cell(I)%corner(1) + cell_width/2.D0)**2 + &
!              (cell(I)%corner(2) + cell_width/2.D0)**2 + &
!              (cell(I)%corner(3) + cell_width/2.D0)**2)
!   print*, r, model_grid(cell(I)%model_index)%rwind, R_inf,  model_grid(cell(I)%model_index)%rho
! END DO
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
  print*, 'we will read a model from Jiri Krticka program...'
  OPEN(UNIT=12,status='old',FILE='jikrmodel.dat')
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
   IF (I == maxrows) THEN
    print*, 'Subroutine read_1D_model:'
    print*, 'Error: Maximum number of records exceeded...'
    print*, 'Exiting program now...'
    STOP
   END IF
   n_modelgrid = n_modelgrid + 1
  END DO
   ALLOCATE (model_grid(n_modelgrid + 1))
   REWIND(11)
  DO I=1,n_modelgrid
   READ(11,*) indexg, r, velo, dens, temp, junk, junk
     model_grid(I)%rwind = r  !  * R_star
     model_grid(I)%vel = velo * 1.E2
     model_grid(I)%rho = dens
     model_grid(I)%T = temp ! should be temp 
     model_grid(I)%J = 0.D0 
     model_grid(I)%assoc_cells = 0
!     print*, 'testing model grid...'
!     print*, model_grid(I)%rwind, model_grid(I)%vel, &
!        model_grid(I)%rho, model_grid(I)%T, model_grid(I)%J, &
!        model_grid(I)%assoc_cells
     !IF (I /= 1) print*, 'delta r: ', model_grid(I)%rwind - model_grid(I-1)%rwind
     ALLOCATE (model_grid(I)%grid_comp(n_elements))
     DO J = 1, n_elements      
        numbions = elements(J)%nions
        ALLOCATE (model_grid(I)%grid_comp(J)%grid_ion(numbions))
        atom_number = elements(J)%atom_number
        model_grid(I)%grid_comp(J)%abund = 1.D0
        !Calculate total number density for included species
        !tot_nd = model_grid(I)%grid_comp(J)%abund / elements(J)%atom_mass 
        !model_grid(I)%grid_comp(J)%numb_den = tot_nd
     END DO
   END DO
  CLOSE(11)
  R_inf  = model_grid(n_modelgrid)%rwind
  V_inf  = model_grid(n_modelgrid)%vel
  ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
  ! All cells out of model grid set to 0 and associate to n_modelgrid. 
  ! Other cells will obtainde particular values with memory
  model_grid(n_modelgrid+1)%rwind = 0.D0
  model_grid(n_modelgrid+1)%vel   = 0.D0
  model_grid(n_modelgrid+1)%rho   = 0.D0     
 CASE DEFAULT
  print*, 'the choice of the variable inputModel = ', inputModel, 'is not known...'
  STOP 'ENDING PROGRAM NOW...'
 END SELECT
  END SUBROUTINE read_1D_model
