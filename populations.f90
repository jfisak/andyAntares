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
 
 pop_number = graund_level_pop * g_stat / g_gstat * &
        exp(e_exc / BOLK / model_grid(model_cell)%T ) / &
        model_grid(model_cell)%volume
! print*, 'populations: pop_number = ', pop_number
! NLTE approximation
CASE(1)
 STOP 'NLTE is not supported yet'
CASE DEFAULT
 STOP 'populations: this choice is not known'
END SELECT
END SUBROUTINE populations
