! subroutine for reading a 2D model grid
SUBROUTINE read_2D_model()

USE types

IMPLICIT NONE
 
 ! variables for reading from a file
 INTEGER                                :: maxrows, ios
 CHARACTER(20)                          :: junk
 ! loop variables
 INTEGER                                :: I, J
 ! readen physical quantities
 DOUBLE PRECISION                       :: radius, perpend, dens, velrad, velang, temp
 INTEGER                                :: atom_number, numbions

 SELECT CASE (inputModel)
  ! reading data from Petr Kurfurst disc model
  ! these files are in this form:
  ! 1. radius / m
  ! 2. coordinate in perpendicular direction to the
  !     disc plane / m
  ! 3. density / kg/m^3
  ! 4. radial velocity / m/s
  ! 5. angular velocity / m/s
  ! 6. temperature / K
  CASE(1)
   print*, 'we will read input input data from Petr Kurfurst model of stellar disc'
   ! firstly we calculate number of rows in the file
   n_modelgrid = 0
   OPEN(UNIT=13,status='old', FILE='disc_model.dat')
    DO I = 1, maxrows
     READ(13,*,IOSTAT = ios) junk
     if(ios /= 0) EXIT
     if(I == maxrows) THEN
      print*, 'maximum number of records exceeded in subroutine read_2d_model'
      print*, 'exiting program now...'
      STOP
     end if
    n_modelgrid = n_modelgrid + 1
   END DO
   ALLOCATE ( model_grid(n_modelgrid + 1))
   REWIND(13)
   DO I = 1, n_modelgrid
    READ(13,*) radius, perpend, dens, velrad, velang, temp
    model_grid(I)%rwind = radius
    model_grid(I)%zwind = perpend
    model_grid(I)%vel = velrad
    model_grid(I)%velang = velang
    model_grid(I)%T = temp
    model_grid(I)%J = 0.D0
    model_grid(I)%assoc_cells = 0
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
   CLOSE(13)
   ! calculation of the stellar radius and Rinf
   R_star = MINVAL(model_grid(:)%rwind)
   R_inf = MAXVAL(model_grid(:)%rwind + model_grid(:)%rwind/1.D2)
  CASE DEFAULT
   STOP 'unknown type of 2D model'
 END SELECT

END SUBROUTINE read_2D_model
