! Read 1D model data and allocate that data to the corresponding values of teh modGrid cells
! the model is 
!
! INPUT: NONE
! OUTPUT: NONE
SUBROUTINE read_1D_average2D() 


USE types
USE constants

IMPLICIT NONE    

INTEGER                                   :: ind_I, ind_J, numbions, atom_number
! for reading from files
DOUBLE PRECISION                          :: junk
INTEGER, PARAMETER                        :: maxrows = 6000000
DOUBLE PRECISION                          :: radius, velo_r, velo_theta, dens, temp
! DOUBLE PRECISION, DIMENSION(n_elements)   :: massfrac
CHARACTER(filename_length)                             :: modelfile
! variables which are not needed in the code
!DOUBLE PRECISION                          :: delta_r, delta, delta2, tot_nd, tot_md
! (2) PoWR model
DOUBLE PRECISION, PARAMETER                :: meanAtMass = 1.33
INTEGER                                         :: reading_grid
INTEGER, DIMENSION(1)                           :: ind_min, ind_max
DOUBLE PRECISION                                :: unit_length, unit_velocity, unit_density
DOUBLE PRECISION                                :: min_radius


modelfile=TRIM(inputmodelFile)
unit_length = 12.64759321736591 * const_Rsun
unit_velocity = 1.D8
unit_density = 1.41314878888978971775872825028550241D-0006

! case
! (1)
! (2) PoWR testing model
! (3) model by Araya

 OPEN (UNIT=11, FILE=modelfile)
 READ(11,*) T_eff
 READ(11,*) R_star
 min_radius = 1.5 * R_star
 ! READ(11,*) R_inf
 ! READ(11,*) V_inf
 ! READ(11,*) V_inf
!  READ(11,*) M_dot
 ! READ(11,*) n_modelgrid
 write(*,*) 'read_1D_average2D: T_eff = ', T_eff
 if(T_eff < 5e3) THEN
  write(*,*) 'T_eff = ', T_eff
  STOP 'effective temperature is too low...'
 end if
 
 
 add_mg = 2
 
 ! Allocate array for model grid structure.
 ! Cell n_modelgrid+1 is associated to propagation grid cells 
 ! which have no counterpart on the modelgrid  
  DO 
   READ(11, *, iostat = reading_grid) radius, junk
   IF(reading_grid /= 0) EXIT
   IF(radius * unit_length < min_radius) cycle
   n_modelgrid = n_modelgrid + 1
  END DO
  ALLOCATE(model_grid(n_modelgrid + add_mg))
  REWIND(11)
  READ(11,*) junk
  READ(11,*) junk
  ! READ(11,*) junk
  ! READ(11,*) junk
  ind_I = 0
  DO 
     ! Maybe better to calculate at the midle of the grid cell rather then at the outer boundary 
     IF(ind_I == n_modelgrid) EXIT
     READ(11,*) radius, velo_r, velo_theta, dens, temp!, massfrac
     IF(radius * unit_length < min_radius) cycle
     ind_I = ind_I + 1
     model_grid(ind_I)%rwind = radius * unit_length
     model_grid(ind_I)%vel = velo_r * unit_velocity
     model_grid(ind_I)%velang = velo_theta * unit_velocity
     model_grid(ind_I)%rho = dens * unit_density
     model_grid(ind_I)%T = temp ! should be temp 
     model_grid(ind_I)%J = 0.D0 
     model_grid(ind_I)%assoc_cells = 0
     ! write(*,'(A23, d14.5)') 'read_1D_model: rwind = ', r * R_star

     ALLOCATE (model_grid(ind_I)%grid_comp(n_elements))
     DO ind_J = 1, n_elements      
        numbions = elements(ind_J)%nions
        ALLOCATE (model_grid(ind_I)%grid_comp(ind_J)%grid_ion(numbions))
        atom_number = elements(ind_J)%atom_number
        model_grid(ind_I)%grid_comp(ind_J)%abund = elements(ind_J)%abundance
     END DO
  END DO
  CLOSE(11)

  R_inf = MAXVAL(model_grid(1:n_modelgrid)%rwind)
  ind_max(1) = MAXLOC(model_grid(:)%rwind,1, MASK=(model_grid(:)%rwind <= R_inf))
  ind_min(1) = MINLOC(model_grid(:)%rwind,1, MASK=(model_grid(:)%rwind >= R_star))
  V_inf = model_grid(ind_max(1))%vel
  V_star = model_grid(ind_min(1))%vel

  write(*,*) 'read_1D_average2D: R_star = ', R_star, ' R_inf = ', R_inf

  ! setting properties
  ! definition of the outerspace index
  outerspace_index = n_modelgrid + 1
  photosphere_index = n_modelgrid + 2


  ! Dummy cell to associate to propagation grid cells which have no representation on the model grid.
  ! All cells out of model grid set to 0 and associate to n_modelgrid. 
  ! Other cells will obtainde particular values with memory
  model_grid(outerspace_index)%rwind = 0.D0
  model_grid(outerspace_index)%vel   = 0.D0
  model_grid(outerspace_index)%rho   = 0.D0     
  model_grid(photosphere_index)%rwind = 0.D0
  model_grid(photosphere_index)%vel   = 0.D0
  model_grid(photosphere_index)%rho   = 0.D0     

! the model cell widths
! now it is a pont in the center of two neighbouring model cells
! the boundaries for the widths calculation
!  ALLOCATE(boundaries(n_modelgrid + 1))
!  boundaries(1) = R_star
!  DO ind_I = 1, n_modelgrid - 1
!   boundaries(ind_I + 1) = (model_grid(ind_I)%rwind + model_grid(ind_I + 1)%rwind)/2.0
!  END DO
!  boundaries(n_modelgrid + 1) = R_inf
!  DO ind_I = 1, n_modelgrid
!   rPrev = boundaries(ind_I)
!   rAct = boundaries(ind_I + 1)
!   width = rAct - rPrev
!   model_grid(ind_I)%width = width
!   ! write(*,*) 'read_1D_model: ind_I = ', I, ' width = ', width
!  END DO
  END SUBROUTINE read_1D_average2D
