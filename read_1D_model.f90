  SUBROUTINE read_1D_model() 

  ! Read 1D model data and allocet that data to the corresponding values of teh model grid cells

  USE types

  IMPLICIT NONE    

  INTEGER                                   :: I, J, M, numbions, indexg, atom_number
  DOUBLE PRECISION                          :: delta_r, delta, delta2, tot_nd, tot_md
  DOUBLE PRECISION                          :: r, velo, dens, temp
  DOUBLE PRECISION, DIMENSION(n_elements)   :: massfrac


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
     model_grid(I)%T = 5000. ! should be temp 
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

  
  END SUBROUTINE read_1D_model
