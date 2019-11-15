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
DOUBLE PRECISION                :: ground_level_pop, g_gstat, g_stat, e_exc
DOUBLE PRECISION                :: abund, rho, atom_mass
DOUBLE PRECISION, PARAMETER     :: minpop = 1.D-50

IF(indexe <= 0) THEN
 write(*,*) 'populations: indexe < 0'
 CALL abort()
END IF
IF(indexi <= 0) THEN
 write(*,*) 'populations: indexi < 0'
 CALL abort()
END IF
IF(level <= 0) THEN
 write(*,*) 'populations: level < 0'
 CALL abort()
END IF

SELECT CASE(nlte)
! LTE approximation
CASE(0)
 ! obvious
 ground_level_pop = model_grid(model_cell)%grid_comp(indexe)%grid_ion(indexi)%gl_pop
 ! statistical weight of the ground state
 g_gstat = elements(indexe)%ions(indexi)%levels(1)%stat_waight
 ! statistical weight of the excited state
 g_stat = elements(indexe)%ions(indexi)%levels(level)%stat_waight
 ! excitation energy
 e_exc = elements(indexe)%ions(indexi)%levels(level)%exci_energy - &
        MINVAL(elements(indexe)%ions(indexi)%levels(:)%exci_energy)
 ! write(*,*) 'populations: e_exc = ', e_exc
 IF(e_exc < 0.D0) STOP 'populations: e_exc < 0'
 rho = model_grid(model_cell)%rho
 abund = model_grid(model_cell)%grid_comp(indexe)%abund
 atom_mass = elements(indexe)%atom_mass
 
 pop_number = 1.E7 * ground_level_pop * g_stat / g_gstat * &
        exp(-e_exc / BOLK / model_grid(model_cell)%T )! * &
 IF(pop_number < minpop) pop_number = 1.D-50
CASE(1)
 STOP 'NLTE is not supported yet'
CASE DEFAULT
 STOP 'populations: this choice is not known'
END SELECT
END SUBROUTINE populations
