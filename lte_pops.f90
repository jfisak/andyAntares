SUBROUTINE lte_pops(indexe, indexi)

USE types
USE constants

IMPLICIT NONE

INTEGER                                 :: indexe, indexi, gridcell
INTEGER, DIMENSION(1)                   :: indexl0
DOUBLE PRECISION                        :: el_nd, temp
DOUBLE PRECISION                        :: frac
DOUBLE PRECISION                        :: gl_pop, N_jk, U

! write(*,*) 'lte_pops: indexe = ', indexe, ' indexi = ', indexi
DO gridcell = 1, n_modelgrid
 IF (model_grid(gridcell)%assoc_cells > 0) THEN
  el_nd = model_grid(gridcell)%e_dens
  temp = model_grid(gridcell)%T
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
  indexl0(1) = MINLOC(elements(indexe)%ions(indexi)%levels(:)%exci_energy,1)
  gl_pop = ( elements(indexe)%ions(indexi)%levels(indexl0(1))%stat_waight * N_jk ) /  U 
  ! IF(indexe == 1 .AND. indexi == 1) write(*,*) 'update_grid: H I = ', frac
  ! IF(indexe == 1 .AND. indexi == 2) write(*,*) 'update_grid: H II = ', frac
  ! write(*,*) 'update_grid: gl_pop = ', gl_pop, ' N_jk = ', N_jk,&
  !  ' U = ', U, ' temp = ', temp
  model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop = gl_pop
  IF(N_jk > 1.D-40) THEN
   model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = N_jk
  ELSE
   model_grid(gridcell)%grid_comp(indexe)%grid_ion(indexi)%tot_pop = 0.E0
  END IF
  ! write(*,*) 'update_grid: indexe = ', indexe, ' indexi = ', indexi, &
  !  ' tot_pop = ', N_jk * frac
  IF(N_jk  > 1.D20) THEN
   write(*,*) 'update_grid: indexe = ', indexe, ' indexi = ', indexi, &
    ' tot_pop = ', N_jk 
   STOP 'lte_pops: suspiciously large number'
  END IF
 END IF
END DO

END SUBROUTINE lte_pops

