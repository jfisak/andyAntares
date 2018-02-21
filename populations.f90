! this subroutine calculate population for the
! given level of the given ion
SUBROUTINE populations(indexe, indexi, level, model_cell, pop_number)
USE types
IMPLICIT NONE

! input variables
INTEGER                         :: indexe, indexi, level, model_cell
! output variables
DOUBLE PRECISION                :: pop_number
! local variables for LTE
DOUBLE PRECISION                :: graund_level_pop, g_gstat, g_stat, e_exc
DOUBLE PRECISION                :: abund, rho, atom_mass
DOUBLE PRECISION, PARAMETER     :: minpop = 1.D-40

SELECT CASE(nlte)
! LTE approximation
CASE(0)
 ! obvious
 graund_level_pop = model_grid(model_cell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop
 ! statistical weight of the ground state
 g_gstat = elements(indexe)%ions(indexi)%levels(1)%stat_waight
 ! statistical weight of the excited state
 g_stat = elements(indexe)%ions(indexi)%levels(level)%stat_waight
 ! excitation energy
 e_exc = elements(indexe)%ions(indexi)%levels(level)%exci_energy - &
        elements(indexe)%ions(indexi)%levels(1)%exci_energy
 rho = model_grid(model_cell)%rho
 abund = model_grid(model_cell)%grid_comp(indexe)%abund
 atom_mass = elements(indexe)%atom_mass
 ! print*, 'populations: indexe = ', indexe, ' indexi = ', indexi, 'populations: e_exc = ', e_exc, &
 ! ' g_stat = ', g_stat, ' g_gstat = ', g_gstat, ' graund_level_pop = ', graund_level_pop
 
 pop_number = graund_level_pop * g_stat / g_gstat * &
        exp(-e_exc / BOLK / model_grid(model_cell)%T )! * &
        ! rho * abund / atom_mass
 IF(pop_number < minpop) pop_number = 1.D-40
! print*, 'populations: rho = ', rho, ' abund = ', abund, ' atom_mass = ', atom_mass
! print*, 'populations: pop_number = ', pop_number
! NLTE approximation
CASE(1)
 STOP 'NLTE is not supported yet'
CASE DEFAULT
 STOP 'populations: this choice is not known'
END SELECT
END SUBROUTINE populations
