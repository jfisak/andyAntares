SUBROUTINE update_grid(iteration)

  ! Calculate electron number density, population number of the ground level
  ! and total population number for every model grid cell for given composition
  ! and corresponding ionization stages.
  USE types

  IMPLICIT NONE    

  INTEGER             :: gridcell, indexe, indexi, numb_ions, iteration
  DOUBLE PRECISION    :: el_nd, temp, frac, U, N_jk, gl_pop
  INTEGER             :: max_n_dcell
  INTEGER             :: I
INTEGER, DIMENSION(1)   :: indexl0

  
  max_n_dcell = SIZE(dyn_cell)
   print*, 'updating grid'
  DO gridcell = 1, n_modelgrid
   !print*, 'update_grid: volume: model cell = ', gridcell, ' volume = ', volume

    IF (model_grid(gridcell)%assoc_cells .GT. 0) THEN
      IF (iteration .EQ. 1) THEN
        ! Calculate electron number density for every model grid cell gridcell
        CALL find_e_nd(gridcell, el_nd)
      ELSE
        ! Energy density contribeted to the model grid cell 
         model_grid(gridcell)%J = model_grid(gridcell)%J / model_grid(gridcell)%volume / (4 * pi)
         temp = (model_grid(gridcell)%J * pi / sigma )**(1./4.) 
!         print*, 'model cell: ', gridcell, ' temperature = ', temp, ' flux = ', model_grid(gridcell)%J, &
!                ' volume = ', model_grid(gridcell)%volume
!         print*, 'update_grid: temperature: I = ', I, ' T = ', temp
         model_grid(gridcell)%T = temp
         ! Calculate electron number density for every model grid cell gridcell
         !print*, gridcell, model_grid(gridcell)%J, temp
         CALL find_e_nd(gridcell, el_nd)
         model_grid(gridcell)%J = 0.D0   
      END IF
      model_grid(gridcell)%e_dens = el_nd 
      temp = model_grid(gridcell)%T

      !     print*, 'temp and e_nd:', gridcell,  model_grid(gridcell)%rho, temp, el_nd/6.1D14
      DO indexe = 1, n_elements
        numb_ions = elements(indexe)%nions
        DO indexi = 1, numb_ions
           ! IF(indexi == numb_ions) THEN
           !  model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop = &
           !   model_grid(gridcell)%grid_comp(indexe)%abund * model_grid(gridcell)%e_dens
           !  model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = &
           !   model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop
           !  CYCLE 
           ! END IF
           ! Calculate fraction (frac) of element indexe in ionization stage indexi
           ! relative to the total number of atoms of this element at given 
           ! electron numb.density and temperature
           CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
           ! Total population number of the element indexe in ionization stage indexi
           ! and particular gridcell (total number of atoms in particular ionization stage)
           N_jk = frac * model_grid(gridcell)%rho * model_grid(gridcell)%grid_comp(indexe)%abund / elements(indexe)%atom_mass
           ! CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
           ! print*, gridcell, indexe, indexi, N_jk/1d10, frac
           ! Calculate partition function (U) of element indexe in ionization stage 
           ! indexi at given temperature temp
           CALL part_fun(indexe, indexi, temp, U)
           ! Ground level population number (number density of the atom at ground level)
           ! write(*,*) 'update_grid: U = ', U
           indexl0(:) = MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
           gl_pop = ( elements(indexe)%ions(indexi)%levels(indexl0(1))%stat_waight * N_jk ) /  U 
           ! write(*,*) 'update_grid: gl_pop = ', gl_pop, ' N_jk = ', N_jk,&
           !  ' U = ', U, ' temp = ', temp
           model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop = gl_pop
           IF(N_jk > 1.D-40) THEN
            model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = N_jk
           ELSE
            model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = 1.D-40
           END IF
           ! write(*,*) 'update_grid: indexe = ', indexe, ' indexi = ', indexi, &
           !  ' tot_pop = ', N_jk * frac
           IF(N_jk  > 1.D20) THEN
            write(*,*) 'update_grid: indexe = ', indexe, ' indexi = ', indexi, &
             ' tot_pop = ', N_jk 
            STOP 'update_grid: suspiciously large number'
           END IF
        END DO
      END DO
    ENDIF
  END DO
  CALL check_pop()
  STOP 'update_grid: testing'
  CLOSE(3)

  
END SUBROUTINE update_grid
