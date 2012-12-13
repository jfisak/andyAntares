SUBROUTINE update_grid()

  ! Calculate electron number density, population number of the ground level
  ! and total population number for every model grin cell for given composition
  ! and corresponding ionization stages.
  USE types

  IMPLICIT NONE    

  INTEGER             :: gridcell, indexe, indexi, numb_ions
  DOUBLE PRECISION    :: el_nd, temp, frac, U, N_jk, gl_pop

  DO gridcell = 1, n_modelgrid
     ! Calculate electron number density for every model grid cell gridcell
     CALL find_e_nd(gridcell, el_nd)
     model_grid(gridcell)%e_dens = el_nd
     temp = model_grid(gridcell)%T
!     print*, 'temp and e_nd:', gridcell,  model_grid(gridcell)%rho, temp, el_nd/6.1D14
     DO indexe = 1, n_elements
        numb_ions = elements(indexe)%nions
        DO indexi = 1, numb_ions
           ! Calculate fraction (frac) of element indexe in ionization stage indexi
           ! relative to the total number of atoms of this element at given 
           ! electron numb.density and temperature
           CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
           ! Total population number of the element indexe in ionization stage indexi
           ! and particular gridcell (total number of atoms in particular ionization stage)
           N_jk = frac * model_grid(gridcell)%rho * model_grid(gridcell)%grid_comp(indexe)%abund / elements(indexe)%atom_mass
           CALL ionization_fraction(indexe, indexi, temp, el_nd, frac)
!           print*, gridcell, indexe, indexi, N_jk/1d10, frac
           ! Calculate partition function (U) of element indexe in ionization stage 
           ! indexi at given temperature temp
           CALL part_fun(indexe, indexi, temp, U)
           ! Ground level population number (number density of the atom at ground level)
           gl_pop = ( elements(indexe)%ions(indexi)%levels(1)%stat_waight * N_jk ) /  U 
           model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop = gl_pop
           model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = N_jk
        END DO
     END DO
!     stop
  END DO
  
END SUBROUTINE update_grid
